#' estimate P site position
#' @description Estimate P site postion from a subset reads.
#' @rdname estimatePsite
#' @param bamfile an BamFile object.
#' @param CDS output of \link{prepareCDS}
#' @param genome a BSgenome object.
#' @return A list with strart codon position, stop codon position from 5 end of the reads,
#' and summary of positions.
#' @references
#' 1: Bazzini AA, Johnstone TG, Christiano R, Mackowiak SD, Obermayer B, Fleming ES,
#' Vejnar CE, Lee MT, Rajewsky N, Walther TC, Giraldez AJ. Identification of small
#' ORFs in vertebrates using ribosome footprinting and evolutionary conservation.
#' EMBO J. 2014 May 2;33(9):981-93. doi: 10.1002/embj.201488411. Epub 2014 Apr 4.
#' PubMed PMID: 24705786; PubMed Central PMCID: PMC4193932.
#' @import GenomicRanges
#' @importFrom BSgenome getSeq
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments qwidth njunc
#' @importFrom Biostrings vmatchPattern
#' @importFrom BiocGenerics table
#' @importFrom methods as is
#' @importFrom XVector subseq
#' @importFrom IRanges subsetByOverlaps
#' @importClassesFrom Rsamtools BamFile
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "tophat2.danRer10.RPF.chr1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' CDS <- readRDS(system.file("extdata", "sampleCDS.rds",
#'                package="ribosomeProfilingQC"))
#' psite <- estimatePsite(bamfile, CDS, Drerio)
#' bestPsite(psite)

estimatePsite <- function(bamfile, CDS, genome){
  stopifnot(is(bamfile, "BamFile"))
  stopifnot(is(CDS, "GRanges"))
  if(length(CDS$internalPos)!=length(CDS) ||
     length(CDS$isLastExonInCDS)!=length(CDS) ||
     length(CDS$tx_name)!=length(CDS) ||
     length(CDS$gene_id)!=length(CDS)){
    stop("CDS must be output of prepareCDS")
  }
  stopifnot(is(genome, "BSgenome"))
  param <-
    ScanBamParam(what=c("qwidth"),
                 tag=character(0),
                 flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                  isUnmappedQuery=FALSE,
                                  isNotPassingQualityControls = FALSE,
                                  isSupplementaryAlignment = FALSE))

  open(bamfile)
  reads <- readGAlignments(bamfile, param = param)
  close(bamfile)

  reads <- reads[qwidth(reads) %in% c(25:30)]
  reads <- reads[njunc(reads)==0]
  x <- as(reads, "GRanges")
  if(length(intersect(seqlevelsStyle(x), seqlevels(CDS)))==0){
    seqlevelsStyle(x) <- seqlevelsStyle(CDS)[1]
  }
  startpos <- getCondonPosition(x, CDS, genome, TRUE)
  stoppos <- getCondonPosition(x, CDS, genome, FALSE)
  stoppos1 <- c(stoppos[-(1:3)], 0, 0, 0)
  pos <- startpos + stoppos1
  x <- assignReadingFrame(promoters(x, upstream = 0, downstream = 1), CDS)
  x <- table(x$readingFrame)
  x <- names(x)[which.max(x)]
  x <- c("0"=13, "1"=14, "2"=12)[x]
  return(list("start codon position"=startpos,
              "stop codon position"=stoppos,
              summary=pos,
              Psite=x))
}

getCodon <- function(CDS.sub, genome, start=TRUE){
  stopifnot(is(CDS.sub, "GRanges"))
  stopifnot(is(genome, "BSgenome"))
  if(start) {
    end(CDS.sub[strand(CDS.sub)=="+"]) <- start(CDS.sub[strand(CDS.sub)=="+"]) + 2
    start(CDS.sub[strand(CDS.sub)=="-"]) <- end(CDS.sub[strand(CDS.sub)=="-"]) - 2
  }else{
    start(CDS.sub[strand(CDS.sub)=="+"]) <- end(CDS.sub[strand(CDS.sub)=="+"]) - 2
    end(CDS.sub[strand(CDS.sub)=="-"]) <- start(CDS.sub[strand(CDS.sub)=="-"]) + 2
  }
  seq <- getSeq(genome, CDS.sub)
  seq <- table(seq)
  seq <- seq[!grepl("N", names(seq))]
  seq <- names(sort(seq, decreasing = TRUE)[1])
  seq
}

#' @rdname estimatePsite
#' @aliases bestPsite
#' @param psite output of estimatePsite function.
#' @return The Psite postion.
#' @importFrom graphics plot text segments
#' @export
bestPsite <- function(psite){
  op <- par(mfrow=c(3, 1))
  on.exit(par(op))
  best.id <- psite$Psite
  psite <- psite[-4]
  mapply(psite, names(psite), FUN=function(.ele, .name){
    xlim <- if(.name == "stop codon position") c(1, 26) else c(-2, 23)
    plot(.ele, main=.name, ylab="count", xlim = xlim)
  })
  x <- seq(from=best.id-3*10, to=23, by=3)
  x <- x[x>0]
  ymax <- max(psite$summary)
  segments(x0 = x-.5, x1 = x - .5,
           y0 = ymax*.02, y1 = ymax*.98,
           lty = 3, col = "gray30")
  text(x=best.id+.25, y= psite$summary[as.character(best.id)],
       labels = "*", col = "red", cex = 2)
  return(best.id)
}


getCondonPosition <- function(x, CDS, genome, start=TRUE){
  stopifnot(is(x, "GRanges"))
  x <- x[width(x)>24]
  x <- unique(x)
  if(length(x)<1000){
    stop("The length of input reads is too less.")
  }
  stopifnot(is(CDS, "GRanges"))
  stopifnot(is(genome, "BSgenome"))
  if(start) {
    stopifnot(length(CDS$isFirstExonInCDS)==length(CDS))
    CDS.sub <- CDS[CDS$isFirstExonInCDS]
  }else{
    stopifnot(length(CDS$isLastExonInCDS)==length(CDS))
    CDS.sub <- CDS[CDS$isLastExonInCDS]
  }
  x <- subsetByOverlaps(x, CDS.sub)
  seq <- getSeq(genome, x)
  seq <- subseq(seq, start = 1, width = 25)
  mp <- vmatchPattern(getCodon(CDS.sub, genome, start), seq)
  mp <- unlist(mp)
  mp <- table(start(mp))
  idx <- as.character(seq.int(23))
  mp <- mp[match(idx, names(mp))]
  mp[is.na(mp)] <- 0
  names(mp) <- idx
  mp
}

