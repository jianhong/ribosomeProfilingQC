#' Plot reads distribution in genomic elements
#' @description Plot the percentage of reads in CDS, 5'UTR, 3'UTR, introns, and
#'  other elements.
#' @param reads Output of \link{getPsiteCoordinates}
#' @param txdb A TxDb object
#' @param upstreamRegion,downstreamRegion The range for promoter region and
#'  downstream region.
#' @param plot Plot the distribution or not
#' @param precedence If no precedence specified, double count will be enabled,
#' which means that if the reads overlap with both CDS and 5'UTR, both
#' CDS and 5'UTR will be incremented. If a precedence order is specified,
#' for example, if promoter is specified before 5'UTR, then only promoter will
#' be incremented for the same example.  The values could be any combinations
#' of "CDS", "UTR5", "UTR3", "OtherExon", "Intron", "upstream", "downstreama"
#' and "InterGenic", Default=NULL
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @param ... Not use.
#' @return The reads with distribution assignment
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
#' pc.sub <- readsDistribution(pc.sub, txdb, las=2,
#'               precedence=c(
#'               "CDS", "UTR5", "UTR3", "OtherExon",
#'               "Intron", "upstream", "downstream",
#'               "InterGenic"
#'               ))

readsDistribution <- function(reads, txdb,
                              upstreamRegion=3000,
                              downstreamRegion=3000,
                              plot=TRUE,
                              precedence=NULL,
                              ignore.seqlevelsStyle=FALSE,
                              ...){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(reads, "GRanges"))
  if(!is.null(precedence)){
    stopifnot(all(precedence %in% c("CDS", "UTR5", "UTR3", "OtherExon",
                                    "Intron", "upstream", "downstream",
                                    "InterGenic")))
  }
  
  reads <- fixSeqlevelsStyle(reads, txdb, ignore.seqlevelsStyle)
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
  type0 <- type <- cbind(type, upstream, downstream, InterGenic)
  if(!is.null(precedence) && all(precedence %in% colnames(type))){
    cn0 <- colnames(type)
    type0 <- type0[, c(precedence, cn0[!cn0 %in% precedence])]
    for(i in seq_along(cn0)[-1]){
      type0[, i] <- type0[, i] & rowSums(type[, seq.int(i-1), drop=FALSE]) == 0
    }
  }
  per <- 100*colSums(type0)/length(reads)
  if(plot) {
    ggBar(per, ylab="percentage (%)", postfix = "%", xlab="")
  }
  mcols(reads) <- cbind(mcols(reads), type)
  return(reads)
}
