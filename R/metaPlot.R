#' metagene analysis plot
#' @description Plot the average coverage of UTR5, CDS and UTR3.
#' @param UTR5coverage,CDScoverage,UTR3coverage translational efficiency of CDS and UTR3 region.
#' Output of \link{coverageDepth}
#' @param sample character(1). Sample name to plot.
#' @param xaxis what to plot for x-axis.
#' @param bins bins for UTR5, CDS and UTR3.
#' @param ... parameter pass to plot.
#' @return a list contain the data for plot.
#' @importFrom IRanges viewMeans Views
#' @export
#' @examples
#'
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cvgs <- coverageDepth(RPFs[1], RNAs[1], gtf)
#' cvgs.utr3 <- coverageDepth(RPFs[1], RNAs[1], gtf, region="utr3")
#' cvgs.utr5 <- coverageDepth(RPFs[1], RNAs[1], gtf, region="utr5")
#' metaPlot(cvgs.utr5, cvgs, cvgs.utr3, sample=1)
#' }
metaPlot <- function(UTR5coverage, CDScoverage, UTR3coverage, sample,
                     xaxis=c("RPFs", "mRNA"), bins=c(UTR5=100, CDS=500, UTR3=100),
                     ...){
  xaxis <- match.arg(xaxis)
  if(!is.list(UTR5coverage) || !is.list(UTR3coverage) || !is.list(CDScoverage)){
    stop("UTR5coverage, CDScoverage and UTR3coverage must be output of coverageDepth")
  }
  if(!xaxis %in% names(UTR5coverage) ||
     !xaxis %in% names(UTR3coverage) ||
     !xaxis %in% names(CDScoverage)){
    stop("UTR5coverage, CDScoverage and UTR3coverage must be output of coverageDepth.",
         "And RPFs or mRNA must be available.")
  }
  if(!is.numeric(sample)){
    sample <- which(names(CDScoverage[[xaxis]]) %in% sample)
  }
  if(length(sample)>1){
    sample <- sample[1]
    message("Only first sample will be plotted.")
  }
  if(!all(c("UTR5", "CDS", "UTR3") %in% names(bins))){
    stop("bins should be a integer vector with names UTR5, CDS and UTR3.")
  }
  regions <- list(UTR5=UTR5coverage[[xaxis]][[sample]],
                  CDS=CDScoverage[[xaxis]][[sample]],
                  UTR3=UTR3coverage[[xaxis]][[sample]])
  bins <- bins[names(regions)]
  regions <- mapply(regions, bins, FUN=function(.ele, .len){
    .ele[lengths(.ele)>=.len]
  }, SIMPLIFY = FALSE)
  ids <- Reduce(intersect, lapply(regions, names))
  regions <- lapply(regions, `[`, i=ids)
  metagene <- mapply(regions, bins, FUN=function(.ele, .len){
    l <- lengths(.ele)
    ir <- IRanges(1, width = l)
    tile <- tile(ir, n = .len)
    names(tile) <- names(.ele)
    vws <- Views(.ele, tile)
    vms <- viewMeans(vws)
    vms <- do.call(rbind, as.list(vms))
    vms <- colMeans(vms)
  }, SIMPLIFY = FALSE)
  v <- unlist(metagene, use.names = TRUE)
  plot(v, type="l", ylab="mean of coverage", xlab="", xaxt="n", ...)
  at <- cumsum(bins)[-3]
  abline(v=at, lty=3, col="gray30")
  labels <- c("START", "STOP")
  axis(1, at = at, labels = labels, las=3)
  at <- cumsum(c(0, bins))[-4]+bins/2
  axis(1, at = at, labels = names(bins), col.ticks = NA)
  return(invisible(metagene))
}
