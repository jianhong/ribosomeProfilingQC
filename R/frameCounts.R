#' extract counts for gene level or transcript level
#' @description Calculate the reads counts or coverage rate for gene level or transcript level.
#' Coverage is determined by measuring the proportion of in-frame CDS positions with >= 1 reads.
#' @param reads output of \link{assignReadingFrame}.
#' @param level transcript or gene level
#' @param frame0only only count for reading frame 0 or not
#' @param coverageRate calculate for coverage or not
#' @return numeric vector with reads counts.
#' @importFrom methods as is
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#' cnts <- frameCounts(pcs)
#' cnts.gene <- frameCounts(pcs, level="gene")
#' cvg <- frameCounts(pcs, coverageRate=TRUE)
frameCounts <- function(reads, level=c("tx", "gene"), frame0only=TRUE, coverageRate=FALSE){
  level <- match.arg(level)
  stopifnot(is(reads, "GRanges"))
  if(length(reads$tx_name)!=length(reads) || length(reads$position)!=length(reads) ||
     length(reads$posToStop)!=length(reads) || length(reads$readingFrame)!=length(reads) ||
     length(reads$gene_id)!=length(reads)){
    stop("reads must be a result of assignReadingFrame")
  }

  if(frame0only){
    reads <- reads[reads$readingFrame %in% 0]
  }
  if(coverageRate){
    reads$CDSwidth <- reads$position + reads$posToStop + 2 # stop codon
    if(frame0only){
      reads$CDSwidth <- reads$CDSwidth/3
    }
    cnt <- split(reads$position, reads$tx_name)
    cnt <- lapply(cnt, unique)
    cnt <- lengths(cnt)
    cnt <- cnt/reads$CDSwidth[match(names(cnt), reads$tx_name)]
    if(level=="gene"){# get mean for each gene
      #cnt <- split(cnt, reads$gene_id[match(names(cnt), reads$tx_name)])
      #cnt <- sapply(cnt, mean)
      stop("gene level coverage rate is not supported by frameCounts!")
    }
  }else{
    if(level=="tx"){
      cnt <- table(reads$tx_name)
    }else{
      cnt <- table(reads$gene_id)
    }
  }
  cnt
}
