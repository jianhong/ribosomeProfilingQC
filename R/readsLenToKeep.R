#' Get reads length to keep by cutoff percentage
#' @description Set the percentage to filter the reads.
#' @param readsLengthDensity Output of \link{summaryReadsLength}
#' @param cutoff Cutoff value.
#' @return Reads length to be kept.
#' @importFrom methods as is
#' @export
#' @examples
#' reads <- GRanges("chr1", ranges=IRanges(seq.int(100), width=1),
#'                  qwidth=sample(25:31, size = 100, replace = TRUE,
#'                                prob = c(.01, .01, .05, .1, .77, .05, .01)))
#' readsLenToKeep(summaryReadsLength(reads, plot=FALSE))
readsLenToKeep <- function(readsLengthDensity, cutoff=0.8){
  stopifnot(is.table(readsLengthDensity))
  readsLengthDensity <- sort(readsLengthDensity, decreasing = TRUE)
  x <- cumsum(readsLengthDensity)
  x <- x[seq.int(which(x>cutoff)[1])]
  x <- range(as.numeric(names(x)))
  seq(from=min(x), to=max(x))
}
