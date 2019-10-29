#' Translational Efficiency
#' @description calculate translational efficiency.
#' @param x output of \link{getFPKM} or \link{normByRUVs}.
#' if window is set, it must be output of \link{coverageDepth}.
#' @param window numeric(1). window size for maximal counts.
#' @param RPFsampleOrder,mRNAsampleOrder sample order of RPFs and mRNAs. The parameters
#' are used to make sure that the order of RPFs and mRNAs in cvgs is corresponding samples.
#' @param pseudocount the number will be add to sum of reads count to avoid X/0.
#' @param log2 Do log2 transform or not.
#' @return a list with RPFs, mRNA levels and TE as a matrix with translational efficiency
#' @importFrom IRanges RleList IRanges IRangesList viewSums
#' @export
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cnts <- countReads(RPFs, RNAs, gtf, level="gene")
#' fpkm <- getFPKM(cnts)
#' te <- translationalEfficiency(fpkm)
#' }
translationalEfficiency <- function(x, window, RPFsampleOrder, mRNAsampleOrder,
                                    pseudocount=1, log2=FALSE){
  if(!is.list(x)){
    stop("x must be output of getFPKM or normByRUVs or coverageDepth.")
  }
  if(!all(c("RPFs", "mRNA") %in% names(x))){
    stop("x must be output of getFPKM or normByRUVs or coverageDepth and must contain RPFs and mRNA.")
  }
  RPFs <- x[["RPFs"]]
  mRNA <- x[["mRNA"]]
  if(missing(window)){
    if(length(dim(RPFs))!=2 | length(dim(mRNA))!=2){
      stop("x must be output of getFPKM or normByRUVs and must contain RPFs and mRNA.")
    }
    if(missing(RPFsampleOrder)) RPFsampleOrder <- seq.int(ncol(x[["RPFs"]]))
    if(missing(mRNAsampleOrder)) mRNAsampleOrder <- seq.int(ncol(x[["mRNA"]]))
    id <- intersect(rownames(RPFs), rownames(mRNA))
    if(length(id)==0){return(NULL)}
    x[["RPFs"]] <- RPFs[id, RPFsampleOrder, drop=FALSE]
    x[["mRNA"]] <- mRNA[id, mRNAsampleOrder, drop=FALSE]
    if(log2){
      x[["TE"]] <- log2(x[["RPFs"]]+pseudocount)- log2(x[["mRNA"]]+pseudocount)
    }else{
      x[["TE"]] <- (x[["RPFs"]]+pseudocount)/(x[["mRNA"]]+pseudocount)
    }
    return(x)
  }else{
    if(!is(RPFs, "list") | !is(mRNA, "list")){
      stop("x must be output of coverageDepth and must contain RPFs and mRNA.")
    }
    if(missing(RPFsampleOrder)) RPFsampleOrder <- seq.int(length(x[["RPFs"]]))
    if(missing(mRNAsampleOrder)) mRNAsampleOrder <- seq.int(length(x[["mRNA"]]))
    RPFs <- RPFs[RPFsampleOrder]
    mRNA <- mRNA[mRNAsampleOrder]
    if(length(mRNA)!=length(RPFs)){
      stop("The length of sample of mRNA is not identical to the length of RPFs.")
    }
    window <- window[1]
    if(round(window)!=window | window < 3){
      stop("window must be a integer no less than 3.")
    }
    ## check all feature_id
    features <- lapply(RPFs, names)
    stopifnot(all(lengths(features)==length(features[[1]])))
    for(i in seq_along(features)){
      stopifnot(all(features[[i]]==features[[1]]))
    }
    features <- features[[1]]

    features2 <- lapply(mRNA, names)
    stopifnot(all(lengths(features2)==length(features2[[1]])))
    for(i in seq_along(features2)){
      stopifnot(all(features2[[i]]==features2[[1]]))
    }
    features2 <- features2[[1]]
    stopifnot(all(features==features2))
    rm(features2)

    cvg <- mapply(RPFs, mRNA, FUN = function(a, b){
      from <- lapply(lengths(a), function(.l){
        seq.int(.l-window+1)
      })
      from[lengths(a)<=window] <- 1
      to <- lapply(from, function(.f){
        .f + window - 1
      })
      to[lengths(a)<=window] <- lengths(a)[lengths(a)<=window]
      ir <- mapply(from, to, FUN = function(f, t){
        IRanges(start = f, end = t)
      })
      ir <- IRangesList(ir)
      vw.a <- Views(a, ir[names(a)])
      vw.b <- Views(b[names(a)], ir[names(a)])
      sum.a <- viewSums(vw.a, na.rm = TRUE)
      sum.b <- viewSums(vw.b, na.rm = TRUE)
      ratios <- mapply(sum.a, sum.b, FUN = function(r, m){
        if(log2){
          log2(r+pseudocount) - log2(m+pseudocount)
        }else{
          (r+pseudocount)/(m+pseudocount)
        }
      }, SIMPLIFY = FALSE)
      ids <- sapply(ratios, which.max)
      ratios <- mapply(ratios, ids, FUN=function(value, key){
        value[key]
      })
      rpf <- mapply(sum.a, ids, FUN=function(value, key){
        value[key]
      })
      mrna <- mapply(sum.b, ids, FUN=function(value, key){
        value[key]
      })
      list(ratios=ratios, RPFs=rpf, mRNA=mrna)
    }, SIMPLIFY = FALSE)
    x[["RPFs"]] <- do.call(cbind, lapply(cvg, `[[`, i="RPFs"))
    x[["mRNA"]] <- do.call(cbind, lapply(cvg, `[[`, i="mRNA"))
    x[["TE"]] <- do.call(cbind, lapply(cvg, `[[`, i="ratios"))
    return(x)
  }
}
