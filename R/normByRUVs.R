#' Normalization by RUVSeq
#' @description Normalization by RUVSeq:RUVs methods
#' @param counts output of \link{countReads}
#' @param RPFgroup,mRNAgroup groups for RPF and mRNA files
#' @param k The number of factor of unwanted variation to be estimated from the data.
#' See \link[RUVSeq:RUVs-methods]{RUVs}
#' @return normalized counts list
#' @importFrom RUVSeq RUVs makeGroups
#' @importFrom EDASeq betweenLaneNormalization newSeqExpressionSet counts normCounts
#' @export
#' @examples
#' \dontrun{##waiting for EDASeq fix the issue.
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' #RPFs <- dir(path, "RPF.*?.[12].bam$", full.names=TRUE)
#' #RNAs <- dir(path, "mRNA.*?.[12].bam$", full.names=TRUE)
#' #gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' #cnts <- countReads(RPFs, RNAs, gtf, level="gene")
#' cnts <- readRDS(file.path(path, "cnts.rds"))
#' gp <- c("KD1", "KD1", "WT", "WT")
#' norm <- normByRUVs(cnts, gp, gp)
#' }
#'

#### Normalization by RUVSeq
#RUVSeq can be used to remove unwanted variation form RNA-Seq data[@risso2014normalization].
#Here we can also use the power of RUVSeq to normalize the count number before we calculate
#translational efficiency.
#
#```{r}
#gp <- c("KD1", "KD1", "KD2", "KD2", "WT", "WT")
#norm <- normByRUVs(cnts, gp)
#```

normByRUVs <- function(counts, RPFgroup, mRNAgroup=RPFgroup, k=1){
  if(!any(c("RPFs", "mRNA") %in% names(counts))){
    stop("counts must be output of coutReads.")
  }
  if(!missing(RPFgroup)){
    if("RPFs" %in% names(counts)){
      RPFs <- counts$RPFs
      if(length(RPFgroup)!=ncol(RPFs)){
        stop("length of RPFgroup is not identical to number of RPF samples")
      }
      counts$RPFs <- normHelper(RPFs, RPFgroup, k)
    }else{
      warning("RPFgroup is set but can not find RPF reads count.")
    }
  }
  if(!missing(mRNAgroup)){
    if("mRNA" %in% names(counts)){
      mRNA <- counts$mRNA
      if(length(mRNAgroup)!=ncol(mRNA)){
        stop("length of mRNAgroup is not identical to number of mRNA samples")
      }
      counts$mRNA <- normHelper(mRNA, mRNAgroup, k)
    }else{
      warning("mRNAgroup is set but can not find mRNA reads count.")
    }
  }
  counts
}

normHelper <- function(cnt, gp, k){
  x <- as.factor(gp)
  set <- newSeqExpressionSet(as.matrix(cnt),
                             phenoData = data.frame(x,
                                                    row.names = colnames(cnt)))
  set1 <- betweenLaneNormalization(set, which = "upper")
  differences <- makeGroups(xs = x)
  nc <- normCounts(set1)
  if(!any(is.na(nc) | is.infinite(nc))){
    set <- set1
  }
  set <- RUVs(x = set, cIdx = rownames(cnt), k=k, scIdx=differences)
  normCounts(set)
}
