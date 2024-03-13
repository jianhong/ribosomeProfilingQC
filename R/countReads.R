#' Extract counts for RPFs and RNAs
#' @description Calculate the reads counts for gene level or transcript level.
#' @param RPFs Bam file names of RPFs.
#' @param RNAs Bam file names of RNAseq.
#' @param gtf GTF file name for annotation.
#' @param level Transcript or gene level.
#' @param bestpsite numeric(1). P site postion.
#' @param readsLen numeric(1). reads length to keep.
#' @param anchor 5end or 3end. Default is 5end.
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @param ... Parameters pass to \link[Rsubread:featureCounts]{featureCounts} 
#' except isGTFAnnotationFile, GTF.attrType, and annot.ext.
#' @return A list with reads counts.
#' @importFrom txdbmaker makeTxDb
#' @importFrom methods as is
#' @importFrom Rsubread featureCounts
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom Rsamtools BamFile
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' RNAs <- dir(path, "mRNA.*?.[12].bam$", full.names = TRUE)
#' cnts <- countReads(RPFs[1], gtf=gtf, level="gene", readsLen=29)
#' #cnts <- countReads(RPFs[1], RNAs[1], gtf=gtf, level="gene", readsLen=29)

countReads <- function(RPFs, RNAs, gtf, level=c("tx", "gene"),
                       bestpsite=13, readsLen=c(28,29), anchor="5end",
                       ignore.seqlevelsStyle=FALSE,
                       ...){
  stopifnot(is.character(gtf))
  level <- match.arg(level)
  anchor <- match.arg(anchor, choices = c("5end", "3end"))
  gtf <- gtf[1]
  counts <- list()
  if(!missing(RPFs)){
    stopifnot(is.character(RPFs))
    stopifnot(is.numeric(readsLen))
    stopifnot(is.numeric(bestpsite))
    suppressWarnings(suppressMessages(txdb <- makeTxDbFromGFF(gtf)))
    counts[["RPFs"]] <- RPFsCounts(files = RPFs,
                                   txdb = txdb, level = level,
                                   bestpsite = bestpsite,
                                   readsLen = readsLen,
                                   anchor = anchor,
                                   ignore.seqlevelsStyle=ignore.seqlevelsStyle)
  }
  if(!missing(RNAs)){
    stopifnot(is.character(RNAs))
    suppressMessages(
      cnts.RNAs <- featureCounts(files = RNAs,
                                 annot.ext = gtf,
                                 GTF.attrType=ifelse(level=="gene",
                                                     "gene_id",
                                                     "transcript_id"),
                                 isGTFAnnotationFile = TRUE,
                                 ...)
    )
    counts[["mRNA"]] <- cnts.RNAs$counts
    counts[["annotation"]] <- cnts.RNAs$annotation
  }else{
    if(missing(RPFs)){
      stop("RNAs or RPFs is required.")
    }
    counts[["annotation"]] <- getFeatureLen(txdb, level)
  }
  if(!missing(RPFs) && !missing(RNAs)){
    ids <- intersect(rownames(counts[["RPFs"]]), rownames(counts[["mRNA"]]))
    counts[["RPFs"]] <- 
      counts[["RPFs"]][match(ids, rownames(counts[["RPFs"]])), , drop=FALSE]
    counts[["mRNA"]] <- 
      counts[["mRNA"]][match(ids, rownames(counts[["mRNA"]])), , drop=FALSE]
    counts[["annotation"]] <- 
      counts[["annotation"]][match(ids, counts[["annotation"]][, 1]), , 
                             drop=FALSE]
  }
  counts
}

RPFsCounts <- function(files, txdb, level, bestpsite,
                       readsLen, anchor,
                       ignore.seqlevelsStyle=FALSE){
  yieldSize <- 10000000
  cnts <- lapply(files, function(f){
    bamfile <- BamFile(file = f, yieldSize = yieldSize)
    pc <- getPsiteCoordinates(bamfile, bestpsite=bestpsite,
                              anchor = anchor)
    pc.sub <- pc[pc$qwidth %in% readsLen]
    pc.sub <- shiftReadsByFrame(pc.sub, txdb,
                                ignore.seqlevelsStyle=ignore.seqlevelsStyle)
    frameCounts(pc.sub, level=level)
  })
  genes <- unique(unlist(lapply(cnts, names)))
  cnts <- lapply(cnts, `[`, i=genes)
  names(cnts) <- basename(files)
  cnts <- do.call(cbind, cnts)
  rownames(cnts) <- genes
  cnts[is.na(cnts)] <- 0
  cnts
}

getFeatureLen <- function(txdb, level=c("gene", "tx")){
  stopifnot(is(txdb, "TxDb"))
  level <- match.arg(level)
  features <- exons(txdb, columns=c("gene_id", "tx_name"))
  features <-
    switch(level,
           gene={
             f <- rep(features, lengths(features$gene_id))
             mcols(f) <-
               DataFrame(feature_id=unlist(features$gene_id))
             f[!is.na(f$feature_id)]
           },
           tx={
             f <- rep(features, lengths(features$tx_name))
             mcols(f) <- DataFrame(unlist(features$tx_name))
             f[!is.na(f$feature_id)]
           })
  f <- as.data.frame(features)
  f <- split(f, f$feature_id)
  GeneID <- names(f)
  foo <- function(.ele, n){
    paste(as.character(.ele[, n]), collapse = ";")
  }
  Chr <- unlist(lapply(f, foo, n="seqnames"))
  Start <- unlist(lapply(f, foo, n="start"))
  End <- unlist(lapply(f, foo, "end"))
  Strand <- unlist(lapply(f, foo, "strand"))
  Length <- unlist(lapply(f, function(.ele) sum(.ele$width)))
  data.frame(GeneID=GeneID, Chr=Chr, Start=Start, End=End,
             Strand=Strand, Length=Length,
             stringsAsFactors = FALSE)
}
