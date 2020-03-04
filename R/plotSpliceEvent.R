#' plot splice event
#' @description Plot the splice event
#' @param se output of \link{spliceEvent}
#' @param tx_name transcript name.
#' @param coverage coverages of feature region with extensions.
#' Output of \link{coverageDepth}
#' @param group1,group2 the sample names of group 1 and group 2
#' @param cutoffFDR cutoff of FDR
#' @param resetIntronWidth logical value. If set to true,
#' reset the region with no read to minimal width.
#' @importFrom ggplot2 geom_segment geom_rect theme element_blank
#' @export
#' @return a ggplot object.
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' coverage <- coverageDepth(RPFs, gtf=gtf, level="gene",
#'                           region="feature with extension")
#' group1 <- c("RPF.KD1.1", "RPF.KD1.2")
#' group2 <- c("RPF.WT.1", "RPF.WT.2")
#' se <- spliceEvent(coverage, group1, group2)
#' plotSpliceEvent(se, se$feature[1], coverage, group1, group2)
#' }

plotSpliceEvent <- function(se, tx_name, coverage, group1, group2,
                            cutoffFDR=0.05, resetIntronWidth=TRUE){
  stopifnot(length(tx_name)==1)
  if(!is(se, "GRanges") || length(se$feature)==0){
    stop("se must be output of spliceEvent")
  }
  if(!is.list(coverage)){
    stop("coverage must be output of coverageDepth")
  }
  if(!"RPFs" %in% names(coverage)){
    stop("coverage must be output of coverageDepth.",
         "And RPFs must be available.")
  }
  if(missing(group1) || missing(group2)){
    stop("group1 or group2 does not have default value")
  }
  cvg <- coverage[["RPFs"]][["coverage"]]
  gr <- coverage[["RPFs"]][["granges"]]
  se.s <- se[se$feature %in% tx_name]
  se.s <- se.s[se.s$FDR<cutoffFDR]
  cvg.s <- lapply(cvg, function(.ele) .ele[[tx_name]])
  cvg.s <- cvg.s[c(group1, group2)]
  gr <- gr[[tx_name]]
  cvg.s <- lapply(cvg.s, as.numeric)
  cvg.s <- do.call(cbind, cvg.s)
  wid <- width(gr)
  id <- rep(seq_along(wid), wid)
  cvg.s <- rowsum(cvg.s, id)
  if(resetIntronWidth){
    cvg.s0 <- rowSums(cvg.s) == 0
    wid[cvg.s0] <- min(wid)
  }
  ol <- findOverlaps(se.s, gr)
  x0 <- cumsum(c(0, wid))
  x <- rep(x0[-length(x0)], ncol(cvg.s))
  xend <- rep(x0[-1], ncol(cvg.s))
  height <- group <- xmax <- xmin <- ymax <- ymin <- NULL
  df <- data.frame(height=log2(as.numeric(cvg.s)+1),
                   sample = rep(colnames(cvg.s), each = nrow(cvg.s)),
                   group = ifelse(rep(colnames(cvg.s),
                                      each = nrow(cvg.s)) %in% group1,
                                  "group1", "group2"),
                   x = x,
                   xend = xend)
  p <- ggplot(data = df, aes(x=x, xend=xend, y=height, yend=height,
                             color = group, fill = sample)) +
    geom_segment() +
    xlab("") + ylab("log2 transformed reads count") +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  if(length(ol)>0) p <- p +
    geom_rect(data = data.frame(xmin = x[subjectHits(ol)],
                                ymin = -.05,
                                xmax = xend[subjectHits(ol)],
                                ymax = max(df$height)*1.05),
              aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax),
              fill = "#33333333", inherit.aes = FALSE)
  p
}
