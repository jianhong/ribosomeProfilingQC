#' Metaplot of P site distribution
#' @description Metaplot of P site distribution in all the CDS aligned by
#' the start codon or stop codon.
#' @param reads output of \link{assignReadingFrame}
#' @param start plot for start codon or stop codon
#' @param anchor xlim position for plot
#' @param col colors for different reading frame.
#' @return invisible height of the barplot.
#' @importFrom methods as is
#' @importFrom graphics barplot axis legend segments
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#' plotDistance2Codon(pcs)
#' plotDistance2Codon(pcs, start=FALSE)
#'
plotDistance2Codon <- function(reads, start=TRUE, anchor=50,
                               col=c("Frame_0" = "#009E73",
                                     "Frame_1" = "#D55E00",
                                     "Frame_2" = "#0072B2")){
  stopifnot(is(reads, "GRanges"))
  if(length(reads$tx_name)!=length(reads) || length(reads$position)!=length(reads) ||
     length(reads$posToStop)!=length(reads) || length(reads$readingFrame)!=length(reads) ||
     length(reads$gene_id)!=length(reads)){
    stop("reads must be a result of assignReadingFrame")
  }
  data <- if(start) reads$position else reads$posToStop
  position <- table(data)
  idx <- seq.int(anchor)
  xlim <- range(idx)
  xlim <- c(0, 1.5*max(xlim))
  if(!start) {
    xlim <- rev(xlim)
    this.col <- col[c(1,3,2)]
  }else{
    this.col <- col
  }
  idx <- c(0, as.character(seq.int(anchor-1)))
  position <- position[match(idx, names(position))]
  position[is.na(position)] <- .1
  names(position) <- idx
  barplot(height=position,
          ylab="Frequence", xlab=paste("P site relative to", ifelse(start, "start", "stop")),
          col=this.col, xlim = xlim,
          border=NA, space=.5, xaxt = "n")
  idx <- as.numeric(idx)
  label <- idx[idx %% 3 == 0][-1]
  axis(1, at = (label+ifelse(start, 0, 1))*1.5 - .5, labels = label)
  axis(1, at = 1, labels = ifelse(start, "START", "STOP"), las = 3)
  if(start){
    idx <- idx[idx %% 3 == 0] * 1.5 + .25
  }else{
    idx <- idx[idx %% 3 == 1] * 1.5 + .25
  }
  segments(x0 = idx, x1 = idx, y0 = 0, y1 = max(position), lty = 3, col = 'gray80')
  legend("topleft", legend = names(col),
         fill = col, border = col, bg = NA, box.col = NA)
  return(invisible(position))
}
