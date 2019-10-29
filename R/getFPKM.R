#' get FPKM values for counts
#' @description Calculate Fragments Per Kilobase of transcript per Million mapped reads (FPKM) for counts.
#' @param counts output of \link{countReads} or \link{normByRUVs}
#' @param gtf GTF file name for annotation.
#' @param level transcript or gene level.
#' @return a list with FPKMs
#' @export
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cnts <- countReads(RPFs, RNAs, gtf, level="gene")
#' fpkm <- getFPKM(cnts)
#' }
getFPKM <- function(counts, gtf, level=c("gene", "tx")){
  if(!is.list(counts)){
    stop("counts must be output of countReads or normByRUVs")
  }
  if(!all(c("RPFs", "mRNA") %in% names(counts))){
    stop("counts must be output of countReads or normByRUVs")
  }
  if(!"annotation" %in% names(counts)){
    annotation <- getFeatureLen(gtf, level)
  }else{
    annotation <- counts$annotation
  }
  FPKM <- function(cnt, annotation){
    total <- colSums(cnt)
    len <- annotation[match(rownames(cnt), annotation$GeneID), "Length"]
    cnt*10e9/total/len
  }
  fpkm <- list()
  if("RPFs" %in% names(counts)){
    fpkm[["RPFs"]] <- FPKM(counts$RPFs, annotation)
  }
  if("mRNA" %in% names(counts)){
    fpkm[["mRNA"]] <- FPKM(counts$mRNA, annotation)
  }
  fpkm
}
