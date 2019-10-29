#' plot density for each reading frame
#' @description plot density for each reading frame
#' @param reads output of \link{assignReadingFrame}
#' @param density plot density or counts
#' @param col colors for reading frames
#' @return reading frame density
#' @importFrom methods as is
#' @importFrom graphics barplot
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#' plotFrameDensity(pcs)
plotFrameDensity <- function(reads, density=TRUE,
                           col=c("Frame_0" = "#009E73",
                                 "Frame_1" = "#D55E00",
                                 "Frame_2" = "#0072B2")){
  stopifnot(is(reads, "GRanges"))
  if(length(reads$readingFrame)!=length(reads)){
    stop("reads must be a result of assignReadingFrame")
  }
  data <- table(reads$readingFrame)
  if(density) data <- data/sum(data)
  barplot(height=data, xlab="Frame",
          ylab=ifelse(density, "Relative Read Density", "Read Count"),
          border=NA, col=col)
  return(data)
}
