#' Codon usage bias
#' @description Calculate the codon usage for the reads in the identified CDSs.
#' And then compared to the reference codon usage.
#' @param RPFs Bam file names of RPFs.
#' @param gtf GTF file name for annotation or a TxDb object.
#' @param genome A BSgenome object.
#' @param bestpsite P site postion.
#' @param readsLen Reads length to keep.
#' @param anchor 5end or 3end. Default is 5end.
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @param summary Return the summary of codon usage bias or full list.
#' @param ... Parameters pass to
#' \link[txdbmaker:makeTxDbFromGFF]{makeTxDbFromGFF}
#' @return A list of data frame of codon count table if summary is TRUE.
#'  Column 'reads' means the counts by raw reads.
#'  Column 'reference' means the counts by sequence extracted from reference by
#'  the coordinates of mapped reads.
#'  PCR duplicated are not removed for the count table.
#'  Otherwise, return the counts (reads/reference) table for each reads.
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings trinucleotideFrequency
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' library(BSgenome.Drerio.UCSC.danRer10)
#' cb <- codonBias(RPFs, gtf=gtf, genome=Drerio)
codonBias <- function(RPFs, gtf, genome,
                      bestpsite=13,
                      readsLen=c(28,29),
                      anchor="5end",
                      ignore.seqlevelsStyle=FALSE,
                      summary=TRUE,
                      ...){
  yieldSize <- 10000000
  stopifnot(is.logical(summary))
  summary <- summary[1]
  stopifnot(is.character(gtf)||is(gtf, "TxDb"))
  stopifnot(is(genome, 'BSgenome'))
  anchor <- match.arg(anchor, choices = c("5end", "3end"))
  stopifnot(is.character(RPFs))
  stopifnot(is.numeric(readsLen))
  stopifnot(is.numeric(bestpsite))
  txdb <- prepareTxDb(gtf, 'gtf', ...)
  cds <- cdsBy(txdb, by = 'tx', use.names = TRUE)
  refSeq <- getSeq(genome, cds)
  refSeq <- vapply(refSeq, paste, FUN.VALUE=character(1L), collapse='')
  res <- lapply(RPFs, function(f){
    bamfile <- BamFile(file = f, yieldSize = yieldSize)
    pc <- getPsiteCoordinates(
      bamfile, bestpsite=bestpsite,
      anchor = anchor,
      param = ScanBamParam(what=c("qwidth", "seq"),
                           tag=character(0),
                           flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                            isUnmappedQuery = FALSE,
                                            isNotPassingQualityControls = FALSE,
                                            isSupplementaryAlignment = FALSE)))
    pc <- pc[pc$qwidth %in% readsLen]
    pc <- shiftReadsByFrame(pc, txdb,
                                ignore.seqlevelsStyle=ignore.seqlevelsStyle)
    pc <- pc[!is.na(pc$tx_name)]
    if(anchor=="5end"){
      pc <- pc[pc$posToStop>0 & pc$position + pc$Psite>0]
    }else{
      pc <- pc[pc$position>0 & pc$posToStop > pc$Psite]
    }
    
    ## position, relative position from start codon
    ## posToStop, relative position from stop codon
    seqStart <- pc$position %% 3
    refStart <- pc$position - seqStart
    k <- refStart<0
    seqStart[k] <- seqStart[k] - refStart[k]
    refStart[k] <- 0
    seqWidth <- pc$qwidth - seqStart
    seqWidth <- seqWidth - (seqWidth %% 3)
    reads <- substr(pc$seq, seqStart+1, seqStart + seqWidth)
    refReads <- substr(refSeq[pc$tx_name],
                       refStart+1, refStart + seqWidth)
    seqCodonUsage <- trinucleotideFrequency(DNAStringSet(reads), step = 3)
    refCodonUsage <- trinucleotideFrequency(DNAStringSet(refReads), step = 3)
    if(summary){
      ## remove the PCR duplicates or not?
      seqCodonUsage <- colSums(seqCodonUsage)
      refCodonUsage <- colSums(refCodonUsage)
      return(data.frame(reads=seqCodonUsage, reference=refCodonUsage))
    }else{
      codonUsage <- paste(seqCodonUsage, refCodonUsage, sep=',')
      dim(codonUsage) <- dim(seqCodonUsage)
      mcols(pc) <- cbind(mcols(pc), codonUsage)
      return(pc)
    }
  })
  names(res) <- basename(RPFs)
  res
}
