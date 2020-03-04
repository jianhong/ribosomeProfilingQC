#' get P site coordinates
#' @description read in all reads from a bamfile and extract P site coordinates.
#' @param bamfile an BamFile object.
#' @param bestpsite P site postion. See \link{estimatePsite}
#' @param anchor 5end or 3end. Default is 5end.
#' @return A GRanges object with qwidth metadata which indicates the width of reads.
#' @import GenomicRanges
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments cigar cigarNarrow cigarQNarrow GAlignments
#' @importFrom methods as is
#' @importFrom S4Vectors metadata<-
#' @importClassesFrom Rsamtools BamFile
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.WT.1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' pc <- getPsiteCoordinates(bamfile, bestpsite=13)
getPsiteCoordinates <- function(bamfile, bestpsite, anchor='5end'){
  stopifnot(is(bamfile, "BamFile"))
  stopifnot(is.numeric(bestpsite))
  anchor <- match.arg(anchor, choices = c("5end", "3end"))
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
    reads <- c(reads, shiftReads(chunk, shift = offset, anchor=anchor))
  }
  close(bamfile)
  reads <- as(reads, "GRanges")
  reads <- promoters(reads, upstream = 0, downstream = 1)
  reads$Psite <- bestpsite
  if(anchor=="3end") reads$Psite <- reads$qwidth - reads$Psite
  metadata(reads) <- list(totalReads=length(reads))
  reads
}


shiftReads <- function(x, shift=12L, anchor="5end"){
  shift <- shift[1]
  stopifnot(round(shift)==shift)
  stopifnot(is(x, "GAlignments"))
  anchor <- match.arg(anchor, choices = c("5end", "3end"))
  if(shift==0){
    return(x)
  }
  strds <- as.character(strand(x)) == "-"
  cigars <- cigar(x)
  cigars <- as.character(cigarNarrow(cigars))
  if(anchor=="5end"){
    cigars <- cigarQNarrow(cigars,
                           start=ifelse(strds, 1, shift+1),
                           end=ifelse(strds, -shift-1, -1))
  }else{
    l <- mcols(x)$qwidth
    shift <- l - shift
    cigars <- cigarQNarrow(cigars,
                           start=ifelse(strds, 1, shift),
                           end=ifelse(strds, -shift, -1))
  }

  x@cigar <- as.character(cigars)
  x@start <- x@start + attributes(cigars)$rshift
  x
}

