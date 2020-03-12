#' Get splicing events
#' @description Get differentical usage of alternative
#' Translation Initiation Sites,
#' alternative Polyadenylation Sites or
#' alternative splicing sites
#' @param coverage Coverages of feature region with extensions.
#' Output of \link{coverageDepth}
#' @param group1,group2 The sample names of group 1 and group 2
#' @importFrom cluster silhouette
#' @importFrom stats fisher.test dist kmeans p.adjust
#' @export
#' @return A GRanges object of splice events.
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' coverage <- coverageDepth(RPFs, gtf=gtf,
#'                   level="gene", region="feature with extension")
#' group1 <- c("RPF.KD1.1", "RPF.KD1.2")
#' group2 <- c("RPF.WT.1", "RPF.WT.2")
#' se <- spliceEvent(coverage, group1, group2)
#' }
spliceEvent <- function(coverage, group1, group2){
  if(!is.list(coverage)){
    stop("coverage must be output of coverageDepth")
  }
  if(!"RPFs" %in% names(coverage)){
    stop("coverage must be output of coverageDepth.",
         "And RPFs must be available.")
  }
  if(missing(group1) || missing(group2)){
    stop("group1 or group2 does not have default value")
  }
  cvg <- coverage[["RPFs"]][["coverage"]]
  gr <- coverage[["RPFs"]][["granges"]]
  if(any(!group1 %in% names(cvg))){
    stop("All the samples should be in coverage.",
         paste(group1[!group1 %in% names(cvg)], collapse = ", "),
         "is not in names of coverage(",
         paste(names(cvg), collapse = ", "),
         ").")
  }
  if(any(!group2 %in% names(cvg))){
    stop("All the samples should be in coverage.",
         paste(group1[!group2 %in% names(cvg)], collapse = ", "),
         "is not in names of coverage(",
         paste(names(cvg), collapse = ", "),
         ").")
  }
  gp1 <- cvg[group1]
  gp2 <- cvg[group2]
  gp1rd <- Reduce(`+`, gp1)
  gp2rd <- Reduce(`+`, gp2)
  rd <- gp1rd + gp2rd
  rds <- vapply(rd, sum, FUN.VALUE = 0)
  keep <- rds>lengths(gr) #reason? average reads number per exon smaller than 1.
  #message("filtered events: ", sum(!keep))
  cvg <- cvg[keep]
  gr <- gr[keep]
  gp1 <- gp1[keep]
  gp2 <- gp2[keep]
  gp1rd <- gp1rd[keep]
  gp2rd <- gp2rd[keep]

  ## get cluster number
  resolution <- 200
  getK <- function(x, y, r){
    x <- log2(as.numeric(x)+1)
    y <- log2(as.numeric(y)+1)
    r1 <- tile(r, width = resolution)
    r1[[1]]$feature_id <- "aTIS"
    start(r1[[1]][1]) <- start(r)[1]
    end(r1[[1]][length(r1[[1]])]) <- end(r)[1]
    r1[[length(r1)]]$feature_id <- "aPAS"
    start(r1[[length(r1)]][1]) <- start(r)[length(r1)]
    end(r1[[length(r1)]][length(r1[[length(r1)]])]) <- end(r)[length(r1)]
    if(length(r)>2){
      r_new <- c(r1[[1]], r[-c(1, length(r))], r1[[length(r1)]])
    }else{
      r_new <- unlist(r1)
    }
    wid <- width(r_new)
    id <- rep(seq_along(wid), wid)
    ## resample
    xs <- rowsum(cbind(x,y), id)
    xs <- 3*xs/wid

    z <- as.numeric(xs)
    n <- length(unique(z))
    if(n==1) return(NULL)
    diss <- dist(xs)
    nc.max <- min(apply(xs, 2, function(.e) length(unique(.e))), 4)
    if(nc.max==1) return(NULL)
    v <- rep(0, nc.max)
    for (i in seq.int(nc.max)[-1]) {
      if(nrow(xs)==i){
        v[i] <- i
      }else{
        clust <- kmeans(xs, i)
        ss <- silhouette(clust$cluster, diss)
        v[i] <- mean(ss[, 3])
      }
    }
    k <- min(which.max(v), n)
    getKm <- function(a){
      if(length(a)==k){
        order(a)
      }else{
        km <- kmeans(a, centers = k)
        oid <- order(km$centers)
        nid <- seq_along(oid)
        names(nid) <- oid
        nid[as.character(km$cluster)]
      }
    }
    km_x <- getKm(xs[, "x"])
    km_y <- getKm(xs[, "y"])
    idx <- which(km_x!=km_y)
    if(length(idx)==0) return(NULL)
    if(length(idx)!=nrow(xs)){
      mu.idx <- xs[-idx, , drop=FALSE]
      mu.idx <- mu.idx[rowSums(mu.idx)!=0, , drop=FALSE]
      if(nrow(mu.idx)>0){
        mu <- colMeans(xs[-idx, , drop=FALSE])
      }else{
        mu <- colMeans(xs)
      }
    }else{
      mu <- colMeans(xs)
    }
    sel <- xs[idx, , drop=FALSE]
    type <- r_new$feature_id[idx]
    xs1 <- xs
    if(any(r_new$feature_id %in% "aTIS")){
      xs1[r_new$feature_id=="aTIS", "x"] <-
        cumsum(rev(xs[r_new$feature_id=="aTIS", "x"]))/
        sum(r_new$feature_id=="aTIS")
      xs1[r_new$feature_id=="aTIS", "y"] <-
        cumsum(rev(xs[r_new$feature_id=="aTIS", "y"]))/
        sum(r_new$feature_id=="aTIS")
    }
    if(any(r_new$feature_id %in% "aPAS")){
      xs1[r_new$feature_id=="aPAS", "x"] <-
        cumsum(xs[r_new$feature_id=="aPAS", "x"])/
        sum(r_new$feature_id=="aPAS")
      xs1[r_new$feature_id=="aPAS", "y"] <-
        cumsum(xs[r_new$feature_id=="aPAS", "y"])/
        sum(r_new$feature_id=="aPAS")
    }
    if(as.character(strand(r_new)[1])!="-") {
      if(any(r_new$feature_id %in% "aTIS"))
        end(r_new[r_new$feature_id=="aTIS"]) <-
          max(end(r_new[r_new$feature_id=="aTIS"]))
      if(any(r_new$feature_id %in% "aPAS"))
        start(r_new[r_new$feature_id=="aPAS"]) <-
          min(start(r_new[r_new$feature_id=="aPAS"]))
    }else{
      if(any(r_new$feature_id %in% "aPAS"))
        end(r_new[r_new$feature_id=="aPAS"]) <-
          max(end(r_new[r_new$feature_id=="aPAS"]))
      if(any(r_new$feature_id %in% "aTIS"))
        start(r_new[r_new$feature_id=="aTIS"]) <-
          min(start(r_new[r_new$feature_id=="aTIS"]))
    }
    wid <- width(r_new)
    pvalue <- vapply(idx, function(i){
      v <- rbind(xs1[i, , drop=FALSE], mu)*wid[i]
      v <- round(2^v)
      tryCatch(fisher.test(v)$p.value,
               error = function(e){
                 1
               })
    }, FUN.VALUE = 0.0)
    out <- r_new[idx]
    mcols(out) <- data.frame(gp1=sel[, 1], gp2=sel[, 2],
                             mu_gp1=mu[1], mu_gp2=mu[2],
                             pvalue, type, stringsAsFactors = FALSE)
    k <- out$type %in% "aTIS"
    if(any(k)){
      k1 <- which.min(out[k]$pvalue)
      k <- which(k)
      k <- k[k!=k1]
      out <- out[-k]
    }
    k <- out$type %in% "aPAS"
    if(any(k)){
      k1 <- which.min(out[k]$pvalue)
      k <- which(k)
      k <- k[k!=k1]
      out <- out[-k]
    }
    out
  }
  kk <- mapply(getK, gp1rd, gp2rd, gr)
  kk <- kk[lengths(kk)>0]
  kk <- unlist(GRangesList(kk))
  kk$FDR <- p.adjust(kk$pvalue, method = "BH")
  kk$feature <- names(kk)
  type <- kk$type
  kk$type <- NULL
  type[!type %in% c("aPAS", "aTIS")] <- "aSE"
  kk$type <- as.factor(type)
  kk[order(kk$FDR)]
}
