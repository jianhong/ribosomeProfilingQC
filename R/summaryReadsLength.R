#' Summary the reads lengths
#' @description Plot the reads length distribution
#' @param reads Output of getPsiteCoordinates
#' @param widthRange The reads range to be plot
#' @param plot Do plot or not
#' @param ... Not use.
#' @return The reads length distribution
#' @importFrom methods as is
#' @importFrom graphics barplot
#' @export
#' @examples
#' reads <- GRanges("chr1", ranges=IRanges(seq.int(100), width=1),
#'                  qwidth=sample(25:31, size = 100, replace = TRUE,
#'                                prob = c(.01, .01, .05, .1, .77, .05, .01)))
#' summaryReadsLength(reads)
summaryReadsLength <- function(reads, widthRange=c(20:35), plot=TRUE, ...){
  stopifnot(is(reads, "GRanges"))
  stopifnot(length(reads$qwidth)==length(reads))
  distribution <- table(reads$qwidth)/sum(table(reads$qwidth))
  per <- distribution[names(distribution) %in% as.character(widthRange)]*100
  ggBar(per, ylab="percentage (%)", xlab="reads length", draw = plot)
  return(sort(distribution, decreasing = TRUE))
}
