#' Calculate ORFscore
#' @description  To calculate the ORFscore, reads were counnted at each position
#' within the ORF.
#' \deqn{ORFscore = log_2((\sum_{n=1}^{3}\frac{(F_i-\bar{F})^2}{\bar{F}}) + 1)}
#' where \eqn{F_n} is the number of reads in reading frame n,
#' \eqn{\bar{F}} is the total number of reads across all three frames
#' divided by 3.
#' If \eqn{F_1} is smaller than \eqn{F_2} or \eqn{F_3},
#' \eqn{ORFscore = -1 X ORFscore}.
#' @references
#' 1: Bazzini AA, Johnstone TG, Christiano R, Mackowiak SD, Obermayer B,
#' Fleming ES, Vejnar CE, Lee MT, Rajewsky N, Walther TC, Giraldez AJ.
#' Identification of small ORFs in vertebrates using ribosome footprinting and
#' evolutionary conservation.
#' EMBO J. 2014 May 2;33(9):981-93.
#' doi: 10.1002/embj.201488411. Epub 2014 Apr 4.
#' PubMed PMID: 24705786; PubMed Central PMCID: PMC4193932.
#' @param reads Output of \link{getPsiteCoordinates}
#' @return A numeric vector with ORFscore.
#' @importFrom methods as is
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#' ORFscore <- getORFscore(pcs)
getORFscore <- function(reads){
  stopifnot(is(reads, "GRanges"))
  if(length(reads$tx_name)!=length(reads)){
    stop("reads must be a result of getReadingFrame")
  }
  Fs <- split(reads$readingFrame, reads$tx_name)
  Fs <- lapply(Fs, table)
  ORFscore <- lapply(Fs, function(.ele){
    .ele <- .ele[c("0", "1", "2")]
    .ele[is.na(.ele)] <- 0
    names(.ele) <- c("0", "1", "2")
    m <- mean(.ele)
    s <- .ele["0"] < .ele["1"] || .ele["0"] < .ele["2"]
    log2(sum((.ele-m)^2/m)+1) * ifelse(s, -1, 1)
  })
  unlist(ORFscore)
}
