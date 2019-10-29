#' get P site coordinates
#' @description read in all reads from a bamfile and extract P site coordinates.
#' @param bamfile an BamFile object.
#' @param bestpsite P site postion. See \link{estimatePsite}
#' @return A GRanges object with qwidth metadata which indicates the width of reads.
#' @import GenomicRanges
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments cigar cigarNarrow cigarQNarrow GAlignments
#' @importFrom methods as is
#' @importClassesFrom Rsamtools BamFile
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "tophat2.danRer10.RPF.chr1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' pc <- getPsiteCoordinates(bamfile, bestpsite=13)
getPsiteCoordinates <- function(bamfile, bestpsite){
  stopifnot(is(bamfile, "BamFile"))
  stopifnot(is.numeric(bestpsite))
  offset <- round(bestpsite[1]) - 1
  param <-
    ScanBamParam(what=c("qwidth"),
                 tag=character(0),
                 flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                  isUnmappedQuery=FALSE,
                                  isNotPassingQualityControls = FALSE,
                                  isSupplementaryAlignment = FALSE))

  reads <- GAlignments()
  open(bamfile)
  while(length(chunk <- readGAlignments(bamfile, param=param))){
    reads <- c(reads, shiftReads(chunk, shift = offset))
  }
  close(bamfile)
  reads <- as(reads, "GRanges")
  reads <- promoters(reads, upstream = 0, downstream = 1)
  reads$Psite <- bestpsite
  reads
}


shiftReads <- function(x, shift=12L){
  shift <- shift[1]
  stopifnot(round(shift)==shift)
  stopifnot(is(x, "GAlignments"))
  if(shift==0){
    return(x)
  }
  strds <- as.character(strand(x)) == "-"
  cigars <- cigar(x)
  cigars <- as.character(cigarNarrow(cigars))
  cigars <- cigarQNarrow(cigars,
                         start=ifelse(strds, 1, shift+1),
                         end=ifelse(strds, -shift-1, -1))
  x@cigar <- as.character(cigars)
  x@start <- x@start + attributes(cigars)$rshift
  x
}

