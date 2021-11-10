fixSeqlevelsStyle <- function(x, CDS){
  if(length(intersect(seqlevelsStyle(x), seqlevelsStyle(CDS)))==0){
    try_res <- try({seqlevelsStyle(x) <- seqlevelsStyle(CDS)[1]})
    if(inherits(try_res, "try-error")){
      warning(try_res)
      warning("The seqlevelsStyles of bam file and the CDS are different.")
    }
  }
  return(x)
}
