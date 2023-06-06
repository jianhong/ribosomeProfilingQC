#' Plot the distribution of reads in sense and antisense strand
#' @description Plot the distribution of reads in sense and antisense strand to
#' check the mapping is correct.
#' @param reads Output of \link{getPsiteCoordinates}
#' @param CDS Output of \link{prepareCDS}
#' @param col Coloar for sense and antisense strand.
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @param ... Parameter passed to barplot
#' @return A ggplot object.
#' @import GenomicRanges
#' @importFrom methods as is
#' @importFrom graphics barplot
#' @importFrom S4Vectors queryHits
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.WT.1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' pc <- getPsiteCoordinates(bamfile, bestpsite=11)
#' pc.sub <- pc[pc$qwidth %in% c(29, 30)]
#' library(GenomicFeatures)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' txdb <- makeTxDbFromGFF(system.file("extdata",
#'           "Danio_rerio.GRCz10.91.chr1.gtf.gz",
#'           package="ribosomeProfilingQC"),
#'           organism = "Danio rerio",
#'           chrominfo = seqinfo(Drerio)["chr1"],
#'           taxonomyId = 7955)
#' CDS <- prepareCDS(txdb)
#' strandPlot(pc.sub, CDS)
#'
strandPlot <- function(reads, CDS, col=c("#009E73", "#D55E00"),
                       ignore.seqlevelsStyle=FALSE, ...){
  stopifnot(is(reads, "GRanges"))
  stopifnot(is(CDS, "GRanges"))
  reads <- fixSeqlevelsStyle(reads, CDS, ignore.seqlevelsStyle)
  ## reads mapped to sense strand
  ol <- findOverlaps(reads, CDS, ignore.strand=FALSE)
  a <- unique(queryHits(ol))
  ## reads mapped to antisense strand
  reads.rev <- switch.strand(reads)
  ol.anti <- findOverlaps(reads.rev, CDS, ignore.strand=FALSE)
  b <- unique(queryHits(ol.anti))
  b <- b[!b %in% a]
  per <- c(sense=length(a), antisense=length(b))/length(reads)*100
  ggBar(per, ylab="mapping rate (%)", xlab="", fill=col, postfix = "%")
}

switch.strand <- function(x){
  str <- strand(x)
  stopifnot(all(levels(str)==c("+", "-", "*")))
  levels(str) <- c("-", "+", "*")
  strand(x) <- factor(str, levels=c("+", "-", "*"))
  x
}
