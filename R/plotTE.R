#' plot translational efficiency
#' @description Scatterplot of RNA/RPFs level compared to the translational efficiency.
#' @param TE output of \link{translationalEfficiency}
#' @param sample character(1). Sample name to plot.
#' @param xaxis what to plot for x-axis.
#' @param removeZero Remove the 0 values from plots.
#' @param log2 do log2 transform or not.
#' @param breaks.length length of breaks for histogram.
#' @param ... parameters pass to \link[graphics:plot]{plot}.
#' @return invisible data.frame with x, y of points.
#' @importFrom graphics hist barplot plot
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' #RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' #RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' #gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' #cnts <- countReads(RPFs, RNAs, gtf, level="gene")
#' cnts <- readRDS(file.path(path, "cnts.rds"))
#' fpkm <- getFPKM(cnts)
#' te <- translationalEfficiency(fpkm)
#' plotTE(te, 1)

plotTE <- function(TE, sample, xaxis=c("mRNA", "RPFs"),
                   removeZero=TRUE, log2=TRUE, breaks.length=50, ...){
  if(!is.list(TE)){
    stop("TE must be output of translationalEfficiency.")
  }
  if(!any(c("RPFs", "mRNA", "TE") %in% names(TE))){
    stop("TE must be output of translationalEfficiency.")
  }
  if(missing(sample)){
    stop("sample is required.")
  }
  xaxis <- match.arg(xaxis)
  mRNA <- TE$mRNA
  RPFs <- TE$RPFs
  x <- TE[[xaxis]]
  TE <- TE$TE
  if(!is.numeric(sample)){
    sample <- which(colnames(TE) %in% sample)
  }
  if(length(sample)>1){
    sample <- sample[1]
    message("Only first sample will be plotted.")
  }
  x <- x[, sample]
  TE <- TE[, sample]
  RPFs <- RPFs[, sample]
  if(removeZero){
    keep <- RPFs>0 & mRNA>0
    if(sum(keep)<1){
      stop("No data available for plotting.")
    }
    x <- x[keep]
    TE <- TE[keep]
  }
  if(log2){
    x <- log2(x)
    TE <- log2(TE)
  }
  opar <- par(fig=c(0, .75, 0, 1), new=FALSE, mar=c(5.1, 4.1, 4.1, 0))
  on.exit(par(opar))
  dots <- list(...)
  args <- dots
  args$x <- x
  args$y <- TE
  if(length(args$xlab)==0) args$xlab <- paste(xaxis, "level")
  if(length(args$ylab)==0) args$ylab <- "Translational Efficiency"
  do.call(plot, args)
  ylim <- par("usr")[3:4]
  par(fig=c(.75, 1, 0, 1), new=TRUE, mar=c(5.1, 0, 4.1, 2.1))
  yhist <- hist(TE, breaks=seq(ylim[1], ylim[2], length.out = breaks.length),
                plot=FALSE)
  args <- dots
  args$height <- yhist$density
  args$axes <- FALSE
  args$space <- 0
  args$horiz <- TRUE
  args$cex <- NULL
  do.call(barplot, args)
  return(invisible(data.frame(x=x, y=TE)))
}
