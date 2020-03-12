#' Plot reads P site abundance for a specific transcript
#' @description Plot the bundances of P site on a transcript.
#' @param reads Output of \link{assignReadingFrame}
#' @param tx_name Transcript names.
#' @param col Colors for reading frames
#' @return Invisible heights of the barplot.
#' @importFrom methods as is
#' @importFrom graphics barplot legend par
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#'
#' plotTranscript(pcs, c("ENSDART00000152562", "ENSDART00000054987"))
plotTranscript <- function(reads, tx_name,
                           col=c("Frame_0" = "#009E73",
                                 "Frame_1" = "#D55E00",
                                 "Frame_2" = "#0072B2")){
  stopifnot(is(reads, "GRanges"))
  if(length(reads$tx_name)!=length(reads) ||
     length(reads$position)!=length(reads) ||
     length(reads$posToStop)!=length(reads) ||
     length(reads$readingFrame)!=length(reads) ||
     length(reads$gene_id)!=length(reads)){
    stop("reads must be a result of assignReadingFrame")
  }
  l <- length(tx_name)
  op <- par(mfrow = c(ceiling(l/floor(sqrt(l))), floor(sqrt(l))))
  on.exit(par(op))
  heights <- list()
  for(i in tx_name){
    x.sub <- reads[reads$tx_name %in% i]
    if(length(x.sub)<1){
      warning("No reads in ", i)
      heights[[i]] <- numeric()
    }else{
      d <- table(mcols(x.sub)[, c("readingFrame", "position")])
      CDS.size <- x.sub[1]$position + x.sub[1]$posToStop + 3
      at <- as.character(seq.int(CDS.size))
      height <- colSums(d)
      height <- height[at]
      height[is.na(height)] <- 0
      names(height) <- at
      cols <- col[apply(d, 2, function(.ele) which(.ele!=0))]
      names(cols) <- colnames(d)
      cols <- cols[at]
      names(cols) <- at
      barplot(xlim=c(0, CDS.size),
              height=height,
              xlab = paste("Base position relative to", i,"CDS"),
              ylab = "Number of reads",
              col = cols, border=cols,  main = i)
      legend("topleft", legend = names(col),
             fill = col, border = col, bg = NA, box.col = NA)
      heights[[i]] <- height
    }
  }
  return(invisible(heights))
}
