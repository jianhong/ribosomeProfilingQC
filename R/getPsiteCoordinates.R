#' Get P site coordinates
#' @description Extract P site coordinates from a bam file to a GRanges object.
#' @param bamfile A BamFile object.
#' @param bestpsite P site postion. See \link{estimatePsite}
#' @param anchor 5end or 3end. Default is 5end.
#' @param param A ScanBamParam object. Please note the 'qwidth' is required.
#' @return A GRanges object with qwidth metadata which indicates the width
#' of reads.
#' @import GenomicRanges
#' @importFrom Rsamtools ScanBamParam scanBamFlag bamWhat
#' @importFrom GenomicAlignments readGAlignments cigar cigarNarrow
#' cigarQNarrow GAlignments qwidth cigarWidthAlongReferenceSpace
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
getPsiteCoordinates <- function(
    bamfile, bestpsite, anchor='5end',
    param = ScanBamParam(what=c("qwidth"),
                         tag=character(0),
                         flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                          isUnmappedQuery=FALSE,
                                          isNotPassingQualityControls = FALSE,
                                          isSupplementaryAlignment = FALSE))){
  stopifnot(is(bamfile, "BamFile"))
  stopifnot('qwidth' %in% bamWhat((param)))
  stopifnot(is.numeric(bestpsite))
  anchor <- match.arg(anchor, choices = c("5end", "3end"))
  offset <- round(bestpsite[1]) - 1

  reads <- GAlignments()
  open(bamfile)
  on.exit(close(bamfile))
  while(length(chunk <- readGAlignments(bamfile, param=param))){
    chunk <- narrow(chunk) ## remove the soft and Hard clips
    chunk <- chunk[njunc(chunk)==0] ## remove the junctions
    reads <- c(reads, shiftReads(chunk, shift = offset, anchor=anchor))
  }
  close(bamfile)
  on.exit()
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
  x <- x[qwidth(x)>shift & width(x)>shift & 
           cigarWidthAlongReferenceSpace(cigar(x), 
                                         N.regions.removed = TRUE)>shift]
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

