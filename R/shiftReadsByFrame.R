#' shift reads by reading frame
#' @description shift reads P site position by reading frame. After shifting, all reading frame will be set as 0
#' @param reads output of \link{getPsiteCoordinates}
#' @param txdb TxDb object.
#' @return reads with reading frame information
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.chr1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' pc <- getPsiteCoordinates(bamfile, bestpsite=11)
#' pc.sub <- pc[pc$qwidth %in% c(29, 30)]
#' library(AnnotationDbi)
#' txdb <- loadDb(system.file("extdata", "danRer10.chr1.txdb",
#'                package="ribosomeProfilingQC"))
#' pc.sub <- shiftReadsByFrame(pc.sub, txdb)
#'
shiftReadsByFrame <- function(reads, txdb){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(reads, "GRanges"))
  if(length(reads$Psite)!=length(reads)){
    stop("reads must be output of getPsiteCoordinates.")
  }
  if(length(intersect(seqlevelsStyle(reads), seqlevels(txdb)))==0){
    seqlevelsStyle(reads) <- seqlevelsStyle(txdb)[1]
  }
  reads <- assignReadingFrame(reads, txdb=txdb)
  nstr <- as.character(strand(reads)) == "-"
  pstr <- !nstr & !is.na(reads$readingFrame)
  nstr <- nstr & !is.na(reads$readingFrame)
  start(reads[pstr]) <- start(reads[pstr]) - reads[pstr]$readingFrame
  reads[pstr]$position <- reads[pstr]$position - reads[pstr]$readingFrame
  reads[pstr]$posToStop <- reads[pstr]$posToStop + reads[pstr]$readingFrame
  reads[pstr]$Psite <- reads[pstr]$Psite - reads[pstr]$readingFrame
  end(reads[nstr]) <- start(reads[nstr]) + reads[nstr]$readingFrame
  start(reads[nstr]) <- end(reads[nstr])
  reads[nstr]$position <- reads[nstr]$position - reads[nstr]$readingFrame
  reads[nstr]$position <- reads[nstr]$position + reads[nstr]$readingFrame
  reads[nstr]$Psite <- reads[nstr]$Psite - reads[nstr]$readingFrame
  width(reads) <- 1
  reads[!is.na(reads$readingFrame)]$readingFrame <- 0
  reads
}
