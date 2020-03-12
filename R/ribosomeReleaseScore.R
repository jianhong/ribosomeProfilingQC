#' Ribosome Release Score (RRS)
#' @description RRS is calculated as the ratio of translational efficiency in
#'  the CDS with RPFs in the 3'UTR.
#' @param cdsTE,utr3TE Translational efficiency of CDS and UTR3 region.
#'  Output of \link{translationalEfficiency}
#' @param CDSsampleOrder,UTR3sampleOrder Sample order of cdsTE and utr3TE.
#'  The parameters are used to make sure that the order of CDS and UTR3 in TE
#'   is corresponding samples.
#' @param pseudocount The number will be add to sum of reads count to avoid X/0.
#' @param log2 Do log2 transform or not.
#' @return A vector of RRS.
#' @export
#' @examples
#'
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cvgs <- coverageDepth(RPFs, RNAs, gtf)
#' cvgs.utr3 <- coverageDepth(RPFs, RNAs, gtf, region="utr3")
#' TE90 <- translationalEfficiency(cvgs, window = 90)
#' TE90.utr3 <- translationalEfficiency(cvgs.utr3, window = 90)
#' rrs <- ribosomeReleaseScore(TE90, TE90.utr3)
#'}
#'
ribosomeReleaseScore <- function(cdsTE, utr3TE, CDSsampleOrder, UTR3sampleOrder,
                                 pseudocount=0, log2=FALSE){
  if(!is.list(cdsTE) || !is.list(utr3TE)){
    stop("cdsTE and utr3TE must be output of translationalEfficiency")
  }
  if(!all(c("RPFs", "mRNA", "TE") %in% names(cdsTE)) ||
     !all(c("RPFs", "mRNA", "TE") %in% names(utr3TE))){
    stop("cdsTE and utr3TE must be output of translationalEfficiency And RPFs,
         mRNA and TE must be available.")
  }
  cdsTE <- cdsTE$TE
  utr3TE <- utr3TE$TE
  id <- intersect(rownames(cdsTE), rownames(utr3TE))
  if(missing(CDSsampleOrder)) CDSsampleOrder <- seq.int(ncol(cdsTE))
  if(missing(UTR3sampleOrder)) UTR3sampleOrder <- seq.int(ncol(utr3TE))
  stopifnot(length(CDSsampleOrder)==length(UTR3sampleOrder))
  cdsTE <- cdsTE[id, CDSsampleOrder, drop=FALSE]
  utr3TE <- utr3TE[id, UTR3sampleOrder, drop=FALSE]
  if(log2){
    log2(cdsTE+pseudocount)-log2(utr3TE+pseudocount)
  }else{
    (cdsTE+pseudocount)/(utr3TE+pseudocount)
  }
}
