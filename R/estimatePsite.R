#' Estimate P site position
#' @description Estimate P site position from a subset reads.
#' @rdname estimatePsite
#' @param bamfile A BamFile object.
#' @param CDS Output of \link{prepareCDS}
#' @param genome A BSgenome object.
#' @param anchor 5end or 3end. Default is 5end.
#' @param readLen The reads length used to estimate.
#' @return A best P site position.
#' @references
#' 1: Bazzini AA, Johnstone TG, Christiano R, Mackowiak SD, Obermayer B,
#' Fleming ES, Vejnar CE, Lee MT, Rajewsky N, Walther TC, Giraldez AJ.
#' Identification of small ORFs in vertebrates using ribosome footprinting and
#' evolutionary conservation. EMBO J. 2014 May 2;33(9):981-93.
#' doi: 10.1002/embj.201488411. Epub 2014 Apr 4.
#' PubMed PMID: 24705786; PubMed Central PMCID: PMC4193932.
#' @import GenomicRanges
#' @importFrom BSgenome getSeq
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments qwidth njunc narrow
#' @importFrom Biostrings vmatchPattern
#' @importFrom BiocGenerics table
#' @importFrom methods as is
#' @importFrom XVector subseq
#' @importFrom IRanges subsetByOverlaps
#' @importClassesFrom Rsamtools BamFile
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.WT.1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' #library(GenomicFeatures)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' #txdb <- makeTxDbFromGFF(system.file("extdata",
#'  #         "Danio_rerio.GRCz10.91.chr1.gtf.gz",
#'  #         package="ribosomeProfilingQC"),
#'  #         organism = "Danio rerio",
#'  #         chrominfo = seqinfo(Drerio)["chr1"],
#'  #         taxonomyId = 7955)
#' #CDS <- prepareCDS(txdb)
#' CDS <- readRDS(system.file("extdata", "CDS.rds",
#'                            package="ribosomeProfilingQC"))
#' estimatePsite(bamfile, CDS, Drerio)
#'

estimatePsite <- function(bamfile, CDS, genome, anchor='5end',
                          readLen=c(25:30)){
  stopifnot(is(bamfile, "BamFile"))
  anchor <- match.arg(anchor, choices = c("5end", "3end"))
  stopifnot(is(CDS, "GRanges"))
  if(length(CDS$internalPos)!=length(CDS) ||
     length(CDS$isLastExonInCDS)!=length(CDS) ||
     length(CDS$tx_name)!=length(CDS) ||
     length(CDS$gene_id)!=length(CDS)){
    stop("CDS must be output of prepareCDS")
  }
  stopifnot(is(genome, "BSgenome"))
  param <-
    ScanBamParam(what=c("qwidth"),
                 tag=character(0),
                 flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                  isUnmappedQuery=FALSE,
                                  isNotPassingQualityControls = FALSE,
                                  isSupplementaryAlignment = FALSE))

  open(bamfile)
  reads <- readGAlignments(bamfile, param = param)
  close(bamfile)
  
  reads <- reads[njunc(reads)==0] ## remove junctions
  reads <- narrow(reads) ## remove the soft and hard clips
  reads <- reads[qwidth(reads) %in% readLen]
  x <- as(reads, "GRanges")
  x <- fixSeqlevelsStyle(x, CDS)
  if(anchor=="3end"){
    x <- switch.strand(x)
    x <- promoters(x, upstream = 0, downstream = 1)
    x <- switch.strand(x)
  }
  x <- assignReadingFrame(promoters(x, upstream = 0, downstream = 1), CDS)
  x <- table(x$readingFrame)
  x <- names(x)[which.max(x)]
  posMap <- list('5end'=c("0"=13, "1"=12, "2"=14),
                 '3end'=c("0"=-16, "1"=-17, "2"=-15))
  x <- posMap[[anchor]][x]
  as.numeric(x)
}

