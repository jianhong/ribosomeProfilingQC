#' Get FPKM values for counts
#' @description Calculate Fragments Per Kilobase of transcript per
#' Million mapped reads (FPKM) for counts.
#' @param counts Output of \link{countReads} or \link{normByRUVs}
#' @param gtf GTF file name for annotation.
#' @param level Transcript or gene level.
#' @return A list with FPKMs
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' #RPFs <- dir(path, "RPF.*?.[12].bam$", full.names=TRUE)
#' #RNAs <- dir(path, "mRNA.*?.[12].bam$", full.names=TRUE)
#' #gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' #cnts <- countReads(RPFs, RNAs, gtf, level="gene")
#' cnts <- readRDS(file.path(path, "cnts.rds"))
#' fpkm <- getFPKM(cnts)
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
    t(t(cnt)/total)/len*10e9
  }
  fpkm <- list()
  if("RPFs" %in% names(counts)){
    fpkm[["RPFs"]] <- FPKM(counts$RPFs, annotation)
    if(is.null(fpkm$RPFsRawCounts)) fpkm[["RPFsRawCounts"]]=counts$RPFs
  }
  if("mRNA" %in% names(counts)){
    fpkm[["mRNA"]] <- FPKM(counts$mRNA, annotation)
    if(is.null(fpkm$mRNARawCounts)) fpkm[["mRNARawCounts"]]=counts$mRNA
  }
  if(all(c("RPFs", "mRNA") %in% names(fpkm))){
    if(!identical(rownames(fpkm[["RPFs"]]), rownames(fpkm[["RPFs"]]))){
      ids <- intersect(rownames(fpkm[["RPFs"]]), rownames(fpkm[["mRNA"]]))
      fpkm[["RPFs"]] <- 
        fpkm[["RPFs"]][match(ids, rownames(fpkm[["RPFs"]])), , drop=FALSE]
      fpkm[["mRNA"]] <- 
        fpkm[["mRNA"]][match(ids, rownames(fpkm[["mRNA"]])), , drop=FALSE]
      fpkm[["RPFsRawCounts"]] <- 
        fpkm[["RPFsRawCounts"]][match(ids, rownames(fpkm[["RPFsRawCounts"]])), 
                                , drop=FALSE]
      fpkm[["mRNARawCounts"]] <- 
        fpkm[["mRNARawCounts"]][match(ids, rownames(fpkm[["mRNARawCounts"]])),
                                , drop=FALSE]
      fpkm[["RPFs"]][is.na(fpkm[["RPFs"]])] <- 0
      fpkm[["mRNA"]][is.na(fpkm[["mRNA"]])] <- 0
      fpkm[["RPFsRawCounts"]][is.na(fpkm[["RPFsRawCounts"]])] <- 0
      fpkm[["mRNARawCounts"]][is.na(fpkm[["mRNARawCounts"]])] <- 0
      fpkm[["annotation"]] <- 
        fpkm[["annotation"]][match(ids, fpkm[["annotation"]][, 1]), , 
                               drop=FALSE]
    }
  }
  fpkm
}
