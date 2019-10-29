#' summary the reads lengths
#' @description plot the reads length distribution
#' @param reads output of getPsiteCoordinates
#' @param widthRange the reads range to be plot
#' @param plot do plot or not
#' @param ... parameters passed to barplot
#' @return the reads length distribution
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
  barplot(distribution[names(distribution) %in% as.character(widthRange)]*100,
          ylab="percentage (%)", xlab="reads length", ...)
  return(sort(distribution, decreasing = TRUE))
}
