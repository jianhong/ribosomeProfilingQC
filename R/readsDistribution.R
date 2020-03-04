#' plot reads distribution in genomic elements
#' @description Plot the percentage of reads in CDS, 5'UTR, 3'UTR, introns, and
#'  other elements.
#' @param reads output of \link{getPsiteCoordinates}
#' @param txdb an TxDb object
#' @param upstreamRegion,downstreamRegion the range for promoter region and
#'  downstream region.
#' @param plot plot the distribution or not
#' @param ... not use.
#' @return the reads with distribution assignment
#' @importFrom GenomicFeatures cds fiveUTRsByTranscript threeUTRsByTranscript
#' exons genes promoters
#' @import GenomicRanges
#' @importFrom methods as is
#' @importFrom S4Vectors subjectHits queryHits
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
#' pc.sub <- readsDistribution(pc.sub, txdb, las=2)

readsDistribution <- function(reads, txdb,
                              upstreamRegion=3000,
                              downstreamRegion=3000,
                              plot=TRUE, ...){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(reads, "GRanges"))
  if(length(intersect(seqlevelsStyle(reads), seqlevels(txdb)))==0){
    seqlevelsStyle(reads) <- seqlevelsStyle(txdb)[1]
  }
  cdss <- cds(txdb)
  utr5 <- unlist(fiveUTRsByTranscript(txdb))
  utr3 <- unlist(threeUTRsByTranscript(txdb))
  exons <- exons(txdb)
  anno <- GRangesList(CDS=cdss, UTR5=utr5, UTR3=utr3, OtherExon=exons)
  anno.ul <- unlist(anno)
  anno.ul$type <- rep(names(anno), lengths(anno))
  anno.dis <- disjoin(anno.ul, with.revmap=TRUE)
  anno.dis <- anno.dis[lengths(anno.dis$revmap)==1]
  anno.dis$revmap <- NULL
  anno$OtherExon <- anno.dis
  type <- lapply(anno, function(.ele){
    ol <- findOverlaps(reads, .ele, ignore.strand=FALSE)
    seq_along(reads) %in% queryHits(ol)
  })
  genes <- genes(txdb)
  ol <- findOverlaps(reads, genes, ignore.strand=FALSE)
  type <- do.call(cbind, type)
  notInExon <- rowSums(type) == 0
  Intron <- seq_along(reads) %in% queryHits(ol) & notInExon
  Promoters <- promoters(genes, upstream = upstreamRegion, downstream = 0)
  ol <- findOverlaps(reads, Promoters, ignore.strand=FALSE)
  upstream <- seq_along(reads) %in% queryHits(ol) & notInExon
  Downstreams <- promoters(switch.strand(genes),
                           upstream = downstreamRegion,
                           downstream = 0)
  Downstreams <- switch.strand(Downstreams)
  ol <- findOverlaps(reads, Downstreams, ignore.strand=FALSE)
  downstream <- seq_along(reads) %in% queryHits(ol) & notInExon
  type <- cbind(type, Intron)
  InterGenic <- rowSums(type) == 0
  type <- cbind(type, upstream, downstream, InterGenic)
  per <- 100*colSums(type)/length(reads)
  if(plot) {
    ggBar(per, ylab="percentage (%)", postfix = "%", xlab="")
  }
  mcols(reads) <- cbind(mcols(reads), type)
  return(reads)
}
