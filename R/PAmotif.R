#' Metaplot of P site distribution
#' @description Metaplot of P site distribution in all the CDS aligned by
#' the start codon or stop codon.
#' @param reads Output of \link{assignReadingFrame} or \link{shiftReadsByFrame}.
#' @param genome A BSgenome object.
#' @param plot Plot the motif or not.
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @return A \link[motifStack:pcm-class]{pcm} object
#' @importFrom IRanges promoters
#' @importFrom Rsamtools getSeq
#' @importFrom Biostrings consensusMatrix AA_STANDARD translate
#' @importClassesFrom motifStack pcm
#' @importFrom motifStack plotMotifLogo
#' @importFrom grid grid.draw grid.newpage
#' @export
#' @examples
#' pcs <- readRDS(system.file("extdata", "samplePc.rds",
#'                package="ribosomeProfilingQC"))
#' library(BSgenome.Drerio.UCSC.danRer10)
#' #PAmotif(pcs, Drerio)
#'
PAmotif <- function(reads, genome, plot=TRUE, ignore.seqlevelsStyle=FALSE){
  stopifnot(is(reads, "GRanges"))
  if(length(reads$tx_name)!=length(reads) ||
     length(reads$position)!=length(reads) ||
     length(reads$posToStop)!=length(reads) ||
     length(reads$readingFrame)!=length(reads) ||
     length(reads$gene_id)!=length(reads)){
    stop("reads must be a result of assignReadingFrame or shiftReadsByFrame")
  }
  stopifnot(is(genome, "BSgenome"))
  reads <- fixSeqlevelsStyle(reads, genome, ignore.seqlevelsStyle)
  reads <- promoters(reads, upstream = 6, downstream = 12)
  seq <- getSeq(genome, reads)
  seq <- translate(seq)
  mat <- consensusMatrix(seq, baseOnly=TRUE)
  pcm <- new("pcm", mat=mat[AA_STANDARD, ], name="PAsite motif")
  if(plot){
    p <- plotMotifLogo(pcm, draw = FALSE, ic.scale = FALSE, ylab = "percentage")
    child <- p[[1]]$children
    child[[which(grepl("xaxis", names(child)))]]$label <-
      c("-2", "-1", "P",  "A",  "1",  "2")
    child[[which(grepl("xaxis", names(child)))]]$children$labels$label <-
      c("-2", "-1", "P",  "A",  "1",  "2")
    p[[1]]$children <- child
    grid.newpage()
    grid.draw(p)
    return(invisible(pcm))
  }else{
    return(pcm)
  }
}
