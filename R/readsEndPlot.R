#' Plot start/stop windows
#' @description Plot the reads shifted from start/stop position of CDS.
#' @param bamfile A BamFile object.
#' @param CDS Output of \link{prepareCDS}
#' @param toStartCodon What to search: start or end codon
#' @param fiveEnd Search from five or three ends of the reads.
#' @param shift number(1). Search from 5' end or 3' end of given number.
#' if fiveEnd set to false, please set the shift as a negative number.
#' @param window The window of CDS region to plot
#' @param readLen The reads length used to plot
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @return The invisible list with counts numbers and reads in GRanges.
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBamHeader
#' @importFrom GenomicAlignments readGAlignments qwidth njunc
#' @importFrom methods as is
#' @importFrom graphics barplot abline
#' @importFrom IRanges Views viewApply
#' @importFrom S4Vectors metadata `metadata<-`
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.WT.1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' #library(txdbmaker)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' #txdb <- makeTxDbFromGFF(system.file("extdata",
#'  #         "Danio_rerio.GRCz10.91.chr1.gtf.gz",
#'  #         package="ribosomeProfilingQC"),
#'  #         organism = "Danio rerio",
#'  #         chrominfo = seqinfo(Drerio)["chr1"],
#'  #         taxonomyId = 7955)
#' #CDS <- prepareCDS(txdb)
#' CDS <- readRDS(system.file("extdata", "CDS.rds",
#'                            package="ribosomeProfilingQC"))
#' re <- readsEndPlot(bamfile, CDS, toStartCodon=TRUE)
#' readsEndPlot(re$reads, CDS, toStartCodon=TRUE, fiveEnd=FALSE)
#' #re <- readsEndPlot(bamfile, CDS, toStartCodon=FALSE)
#' #readsEndPlot(re$reads, CDS, toStartCodon=FALSE, fiveEnd=FALSE)
#' readsEndPlot(bamfile, CDS, shift=13)
#' #readsEndPlot(bamfile, CDS, fiveEnd=FALSE, shift=-16)
readsEndPlot <- function(bamfile, CDS, toStartCodon=TRUE,
                         fiveEnd=TRUE, shift=0, window=c(-29, 30),
                         readLen=25:30, ignore.seqlevelsStyle=FALSE){
  stopifnot(is(fiveEnd, "logical"))
  stopifnot(is(shift, "numeric"))
  stopifnot(is(window, "numeric"))
  stopifnot(is(readLen, "numeric"))
  stopifnot(is(CDS, "GRanges"))
  if(length(CDS$internalPos)!=length(CDS) ||
     length(CDS$isFirstExonInCDS)!=length(CDS) ||
     length(CDS$isLastExonInCDS)!=length(CDS) ||
     length(CDS$tx_name)!=length(CDS) ||
     length(CDS$gene_id)!=length(CDS)){
    stop("CDS must be output of prepareCDS")
  }
  
  stopifnot(inherits(bamfile, c("BamFile", "GRanges")))
  argg <- as.list(match.call())
  argg$fiveEnd <- NULL
  if(is(bamfile, "GRanges")){
    metadata <- metadata(bamfile)
    which <- metadata$which
    metadata$which <- NULL
    argg$bamfile <- metadata$bamfile
    if(!identical(argg, metadata)){
      stop("The reads are not ready for current parameter set. ",
           "Please try the BamFile again.")
    }
    reads <- bamfile
  }else{
    if(toStartCodon){
      CDS <- CDS[CDS$isFirstExonInCDS]
      which <- promoters(CDS, upstream=abs(window)[1],
                         downstream=abs(window)[2])
    }else{
      CDS <- CDS[CDS$isLastExonInCDS]
      CDS <- switch.strand(CDS)
      which <- promoters(CDS, upstream=abs(window)[1],
                         downstream=abs(window)[2])
      which <- switch.strand(which)
    }
    
    h <- scanBamHeader(bamfile)
    seqs <- h$targets
    which <- fixSeqlevelsStyle(which, names(seqs), ignore.seqlevelsStyle)
    which <- which[as.character(seqnames(which)) %in% names(seqs)]
    seqlevels(which) <-
      seqlevels(which)[seqlevels(which) %in% names(seqs)]
    param <-
      ScanBamParam(what=c("qwidth"),
                   tag=character(0),
                   flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                    isUnmappedQuery=FALSE,
                                    isNotPassingQualityControls = FALSE,
                                    isSupplementaryAlignment = FALSE),
                   which = which
      )
    
    open(bamfile)
    on.exit(close(bamfile))
    reads <- readGAlignments(bamfile, param = param)
    close(bamfile)
    on.exit()
    
    reads <- reads[njunc(reads)==0]
    reads <- narrow(reads) ## remove the soft and hard clips
    reads <- reads[qwidth(reads) %in% readLen]
    reads <- as(reads, "GRanges")
    reads <- fixSeqlevelsStyle(reads, CDS, ignore.seqlevelsStyle)
    which <- fixSeqlevelsStyle(which, CDS, ignore.seqlevelsStyle)
    argg$which <- which
    metadata(reads) <- argg
  }
  
  if(fiveEnd[1]){
    x <- promoters(reads, upstream = 0, downstream = 1)
    if(shift[1] != 0){
      x <- shift(x, shift = shift[1]-1)
    }
  }else{
    x <- switch.strand(reads)
    x <- promoters(x, upstream = 0, downstream = 1)
    x <- switch.strand(x)
    if(shift[1] != 0){
      x <- shift(x, shift = shift[1]+1)
    }
  }
  x <- split(x, strand(x))
  which <- split(which, strand(which))
  heights <- mapply(x, which[names(x)], names(x),
                    FUN=function(.x, .which, .n){
    if(length(.x)==0){
      return(NULL)
    }
    cvg <- coverage(.x)
    w <- split(.which, seqnames(.which))
    cvg.sub <- unlist(lapply(cvg, sum))
    cvg <- cvg[cvg.sub>0]
    seq <- intersect(names(cvg), names(w))
    vws <- Views(cvg[seq], w[seq])
    vws <- lapply(vws, function(.ele) {
      viewApply(.ele[width(.ele)==sum(abs(window)[c(1, 2)])], as.numeric)
    })
    vws <- do.call(cbind, vws)
    at <- seq(-abs(window[1]), abs(window[2]))
    at <- at[at!=0]
    if(.n=='-'){
      at <- rev(at)
    }
    if(length(dim(vws))!=2){
      stop("Not enough data available.")
    }
    height <- rowSums(vws)
    names(height) <- at
    height[order(at)]
  })
  heights <- do.call(cbind, heights)
  heights <- rowSums(heights)
  at <- as.numeric(names(heights))
  barplot(heights, las=3, space = .5, ylab = "counts",
          xlab = paste("distance from", ifelse(fiveEnd, "5'", "3'"),
                       "of reads to",
                       ifelse(toStartCodon, "start", "stop"),
                       "codon"))
  at1 <- which(at==-1)
  abline(v=at1*1.5+.25, lty=4)
  at <- seq((at1 %% 3), to = length(at), by = 3)
  at <- at[at!=at1]
  at <- at * 1.5 + .25
  ymax <- max(heights)
  segments(x0 = at, y0 = 0, x1 = at, y1 = ymax * .9,
           lty = 3, col = "gray80")
  return(invisible(list(reads=reads, heights=heights)))
}

