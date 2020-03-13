#' Start or Stop codon usage
#' @description Calculate the start or stop codon usage for the identified CDSs.
#' @param reads Output of \link{assignReadingFrame}.
#' @param start Calculate for start codon or stop codon.
#' @param genome A BSgenome object.
#' @return Table of codon usage.
#' @importFrom methods as is
#' @importFrom Biostrings getSeq
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#' library(BSgenome.Drerio.UCSC.danRer10)
#' codonUsage(pcs, genome=Drerio)
#' codonUsage(pcs, start=FALSE, genome=Drerio)
codonUsage <- function(reads, start=TRUE, genome){
  stopifnot(is(reads, "GRanges"))
  stopifnot(is.logical(start))
  if(length(reads$tx_name)!=length(reads) ||
     length(reads$position)!=length(reads) ||
     length(reads$posToStop)!=length(reads) ||
     length(reads$readingFrame)!=length(reads) ||
     length(reads$gene_id)!=length(reads)){
    stop("reads must be a result of assignReadingFrame")
  }
  if(start){
    reads <- reads[!is.na(reads$position)]
    reads <- reads[reads$position==0]
  }  else {
    reads <- reads[!is.na(reads$posToStop)]
    reads <- reads[reads$posToStop==0]
  }
  reads <- promoters(reads, upstream = 0, downstream = 3)
  seq <- getSeq(genome, reads)
  sort(table(seq), decreasing = TRUE)
}
