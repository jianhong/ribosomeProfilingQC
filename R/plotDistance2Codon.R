#' Metaplot of P site distribution
#' @description Metaplot of P site distribution in all the CDS aligned by
#' the start codon or stop codon.
#' @param reads Output of \link{assignReadingFrame}.
#' @param start Plot for start codon or stop codon.
#' @param anchor The maximal xlim or (min, max) position for plot.
#' @param col Colors for different reading frame.
#' @return Invisible height of the barplot.
#' @importFrom methods as is
#' @importFrom graphics barplot axis legend segments
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#' plotDistance2Codon(pcs)
#' #plotDistance2Codon(pcs, start=FALSE)
#' #plotDistance2Codon(pcs, anchor=c(-10, 20))
plotDistance2Codon <- function(reads, start=TRUE, anchor=50,
                               col=c("Frame_0" = "#009E73",
                                     "Frame_1" = "#D55E00",
                                     "Frame_2" = "#0072B2")){
  stopifnot(is(reads, "GRanges"))
  stopifnot(is(anchor, "numeric"))
  stopifnot(is.logical(start))
  if(length(reads$tx_name)!=length(reads) ||
     length(reads$position)!=length(reads) ||
     length(reads$posToStop)!=length(reads) ||
     length(reads$readingFrame)!=length(reads) ||
     length(reads$gene_id)!=length(reads)){
    stop("reads must be a result of assignReadingFrame")
  }
  data <- if(start) reads$position else reads$posToStop
  position <- table(data)
  if(length(anchor)==0) {
    stop("Length of anchor must be greater than 0.")
  }
  if(length(anchor)==1){
    anchor <- c(0, anchor)
  }
  idx <- seq(from = anchor[1], to = anchor[2])
  xlim <- range(c(0, length(idx)))
  xlim <- 1.5 * xlim
  if(!start) {
    xlim <- rev(xlim)
    this.col <- col[c(1,3,2)]
  }else{
    this.col <- col
  }
  idx <- as.character(idx)
  position <- position[match(idx, names(position))]
  position[is.na(position)] <- .01
  names(position) <- idx
  barplot(height=position,
          ylab="Frequence",
          xlab=paste("P site relative to", ifelse(start, "start", "stop")),
          col=this.col[abs(as.numeric(names(position)) %% 3) + 1],
          xlim = xlim,
          border=NA, space=.5, xaxt = "n")
  idx <- as.numeric(idx)
  at <- which(idx %% 3 == 0 & idx != 0)
  at <- at
  label <- as.character(idx[at])
  axis(1, at = at*1.5 - .5, labels = label)
  axis(1, at = which(idx==0)*1.5-.5,
       labels = ifelse(start, "START", "STOP"), las = 3)
  if(start){
    at <- which(idx %% 3 == 2) * 1.5 + .25
  }else{
    at <- which(idx %% 3 == 0) * 1.5 + .25
  }
  segments(x0 = at, x1 = at, y0 = 0, y1 = max(position),
           lty = 3, col = 'gray80')
  legend("topleft", legend = names(col),
         fill = col, border = col, bg = NA, box.col = NA)
  return(invisible(position))
}
