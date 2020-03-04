#' filter CDS by size
#' @description filter CDS by CDS size.
#' @param CDS output of preparedCDS
#' @param sizeCutoff cutoff size for CDS.
#' If the size of CDS is less than the cutoff, it will be filtered out.
#' @return an GRanges object with filtered CDS.
#' @export
#' @import GenomicRanges
#' @importFrom methods as is
#' @examples
#' library(GenomicFeatures)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' txdb <- makeTxDbFromGFF(system.file("extdata",
#'           "Danio_rerio.GRCz10.91.chr1.gtf.gz",
#'           package="ribosomeProfilingQC"),
#'           organism = "Danio rerio",
#'           chrominfo = seqinfo(Drerio)["chr1"],
#'           taxonomyId = 7955)
#' CDS <- prepareCDS(txdb)
#' filterCDS(CDS)
filterCDS <- function(CDS, sizeCutoff = 100L){
  sizeCutoff <- sizeCutoff[1]
  stopifnot(is.numeric(sizeCutoff))
  stopifnot(is(CDS, "GRanges"))
  if(length(CDS$internalPos)!=length(CDS) ||
     length(CDS$isLastExonInCDS)!=length(CDS) ||
     length(CDS$tx_name)!=length(CDS) ||
     length(CDS$gene_id)!=length(CDS)){
    stop("CDS must be output of prepareCDS")
  }
  keep <- unique(CDS[CDS$wid.cumsum>=sizeCutoff]$tx_name)
  CDS[CDS$tx_name %in% keep]
}
