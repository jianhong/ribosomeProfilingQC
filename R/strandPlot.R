#' plot the distribution of reads in sense and antisense strand
#' @description plot the distribution of reads in sense and antisense strand to
#' check the mapping is correct.
#' @param reads output of \link{getPsiteCoordinates}
#' @param CDS output of \link{prepareCDS}
#' @param col coloar for sense and antisense strand.
#' @param ... parameter passed to barplot
#' @return ggplot object.
#' @import GenomicRanges
#' @importFrom methods as is
#' @importFrom graphics barplot
#' @importFrom S4Vectors queryHits
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "tophat2.danRer10.RPF.chr1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' pc <- getPsiteCoordinates(bamfile, bestpsite=11)
#' pc.sub <- pc[pc$qwidth %in% c(29, 30)]
#' CDS <- readRDS(system.file("extdata", "sampleCDS.rds",
#'                package="ribosomeProfilingQC"))
#' strandPlot(pc.sub, CDS)
#'
strandPlot <- function(reads, CDS, col=c("#009E73", "#D55E00"), ...){
  stopifnot(is(reads, "GRanges"))
  stopifnot(is(CDS, "GRanges"))
  if(length(intersect(seqlevelsStyle(reads), seqlevels(CDS)))==0){
    seqlevelsStyle(reads) <- seqlevelsStyle(CDS)[1]
  }
  ## reads mapped to sense strand
  ol <- findOverlaps(reads, CDS, ignore.strand=FALSE)
  a <- length(unique(queryHits(ol)))/length(reads)
  ## reads mapped to antisense strand
  reads.rev <- switch.strand(reads)
  ol.anti <- findOverlaps(reads.rev, CDS, ignore.strand=FALSE)
  b <- length(unique(queryHits(ol.anti)))/length(reads)
  per <- c(sense=a, antisense=b)*100
  ggBar(per, ylab="mapping rate (%)", xlab="", fill=col, postfix = "%")
}

switch.strand <- function(x){
  str <- strand(x)
  stopifnot(all(levels(str)==c("+", "-", "*")))
  levels(str) <- c("-", "+", "*")
  strand(x) <- factor(str, levels=c("+", "-", "*"))
  x
}
