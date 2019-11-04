#' plot density for each reading frame
#' @description plot density for each reading frame
#' @param reads output of \link{assignReadingFrame}
#' @param density plot density or counts
#' @param col colors for reading frames
#' @return reading frame density
#' @importFrom methods as is
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
  ggBar(height = data*100, xlab="Frame",
        ylab=ifelse(density, "Relative Read Density (%)", "Read Count"),
        fill = col, postfix = ifelse(density, "%", FALSE))
}
