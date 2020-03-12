#' Plot start/stop windows
#' @description Plot the reads shifted from start/stop position of CDS.
#' @param bamfile A BamFile object.
#' @param CDS Output of \link{prepareCDS}
#' @param toStartCodon What to search: start or end codon
#' @param fiveEnd Search from five or three ends of the reads
#' @param window The window of CDS region to plot
#' @param readLen The reads length used to plot
#' @return The invisible counts numbers.
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBamHeader
#' @importFrom GenomicAlignments readGAlignments qwidth njunc
#' @importFrom methods as is
#' @importFrom graphics barplot abline
#' @importFrom IRanges Views viewApply
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.WT.1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' #library(GenomicFeatures)
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
#' readsEndPlot(bamfile, CDS, toStartCodon=TRUE)
#' readsEndPlot(bamfile, CDS, toStartCodon=TRUE, fiveEnd=FALSE)
#' #readsEndPlot(bamfile, CDS, toStartCodon=FALSE)
#' #readsEndPlot(bamfile, CDS, toStartCodon=FALSE, fiveEnd=FALSE)
readsEndPlot <- function(bamfile, CDS, toStartCodon=TRUE,
                         fiveEnd=TRUE, window=c(-29, 30),
                         readLen=25:30){
  stopifnot(is(bamfile, "BamFile"))
  stopifnot(is(CDS, "GRanges"))
  if(length(CDS$internalPos)!=length(CDS) ||
     length(CDS$isFirstExonInCDS)!=length(CDS) ||
     length(CDS$isLastExonInCDS)!=length(CDS) ||
     length(CDS$tx_name)!=length(CDS) ||
     length(CDS$gene_id)!=length(CDS)){
    stop("CDS must be output of prepareCDS")
  }
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
  if(length(intersect(seqlevelsStyle(names(seqs)),
                      seqlevelsStyle(which)))==0){
    seqlevelsStyle(which) <- seqlevelsStyle(names(seqs))[1]
    which <- which[as.character(seqnames(which)) %in% names(seqs)]
    seqlevels(which) <-
      seqlevels(which)[seqlevels(which) %in% names(seqs)]
  }
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
  reads <- readGAlignments(bamfile, param = param)
  close(bamfile)

  reads <- reads[qwidth(reads) %in% readLen]
  reads <- reads[njunc(reads)==0]
  reads <- as(reads, "GRanges")
  if(length(intersect(seqlevelsStyle(reads), seqlevels(CDS)))==0){
    seqlevelsStyle(reads) <- seqlevelsStyle(CDS)[1]
    seqlevelsStyle(which) <- seqlevelsStyle(CDS)[1]
  }
  if(fiveEnd){
    x <- promoters(reads, upstream = 0, downstream = 1)
  }else{
    x <- switch.strand(reads)
    x <- promoters(x, upstream = 0, downstream = 1)
    x <- switch.strand(x)
  }
  cvg <- coverage(x)
  w <- split(which, seqnames(which))
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
  if(length(dim(vws))!=2){
    stop("Not enough data available.")
  }
  height <- rowSums(vws)
  names(height) <- at
  barplot(height, las=3, space = .5, ylab = "counts",
          xlab = paste("distance from", ifelse(fiveEnd, "5'", "3'"),
                       "of reads to",
                       ifelse(toStartCodon, "start", "stop"),
                       "codon"))
  at1 <- which(at==-1)
  abline(v=at1*1.5+.25, lty=4)
  at <- seq((at1 %% 3), to = length(at), by = 3)
  at <- at[at!=at1]
  at <- at * 1.5 + .25
  ymax <- max(height)
  segments(x0 = at, y0 = 0, x1 = at, y1 = ymax * .9,
           lty = 3, col = "gray80")
  return(invisible(height))
}

