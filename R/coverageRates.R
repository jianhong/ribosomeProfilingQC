#' calculate coverage rate
#' @description calculate coverage rate for RPFs and mRNAs
#' Coverage will be calculated based on best P sites for RPFs and 5'end for RNA-seq.
#' @param cvgs output of \link{coverageDepth}
#' @param RPFsampleOrder,mRNAsampleOrder sample order of RPFs and mRNAs. The parameters
#' are used to make sure that the order of RPFs and mRNAs in cvgs is corresponding samples.
#' @return a list with coverage rate
#' @export
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cvgs <- coverageDepth(RPFs, RNAs, gtf, level="tx")
#' cr <- coverageRates(cvgs)
#' }
coverageRates <- function(cvgs, RPFsampleOrder, mRNAsampleOrder){
  if(!is.list(cvgs)){
    stop("cvgs must be output of coverageDepth.")
  }
  if(!any(c("RPFs", "mRNA") %in% names(cvgs))){
    stop("cvgs must be output of coverageDepth.")
  }
  cr <- lapply(cvgs, function(.ele){
    .ele <- .ele[["coverage"]]
    covered <- lapply(.ele, function(.e){
      sapply(.e>0, sum)/lengths(.e)
    })
    ids <- unique(unlist(lapply(covered, names)))
    covered <- lapply(covered, function(.e){
      .e[ids]
    })
    covered <- do.call(cbind, covered)
    covered[is.na(covered)] <- 0
    rownames(covered) <- ids
    covered
  })

  if("RPFs" %in% names(cr)) cr[["RPFs"]] <- cr[["RPFs"]]*3
  if(all(c("RPFs", "mRNA") %in% names(cr))){
    if(missing(RPFsampleOrder)) RPFsampleOrder <- seq.int(ncol(cr[["RPFs"]]))
    if(missing(mRNAsampleOrder)) mRNAsampleOrder <- seq.int(ncol(cr[["mRNA"]]))
    ids <- intersect(rownames(cr[["RPFs"]]), rownames(cr[["mRNA"]]))
    cr[["coverageRatio"]] <- cr[["RPFs"]][ids, RPFsampleOrder]/cr[["mRNA"]][ids, mRNAsampleOrder]
  }
  cr
}
