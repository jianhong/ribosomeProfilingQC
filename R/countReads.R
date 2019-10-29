#' extract counts for RPFs and RNAs
#' @description Calculate the reads counts for gene level or transcript level.
#' @param RPFs Bam file names of RPFs.
#' @param RNAs Bam file names of RNAseq.
#' @param gtf GTF file name for annotation.
#' @param level transcript or gene level.
#' @param bestpsite numeric(1). P site postion.
#' @param readsLen numeric(1). reads length to keep.
#' @param ... parameters pass to \link[Rsubread:featureCounts]{featureCounts}
#' @return a list with reads counts.
#' @importFrom methods as is
#' @importFrom Rsubread featureCounts
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom Rsamtools BamFile
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cnts <- countReads(RPFs, RNAs, gtf, level="tx")
#' }

countReads <- function(RPFs, RNAs, gtf, level=c("tx", "gene"),
                       bestpsite=13, readsLen=c(28,29),
                       ...){
  stopifnot(is.character(gtf))
  level <- match.arg(level)
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
                                   readsLen = readsLen)
  }
  if(!missing(RNAs)){
    stopifnot(is.character(RNAs))
    suppressMessages(
      cnts.RNAs <- featureCounts(files = RNAs,
                                 annot.ext = gtf,
                                 GTF.attrType=ifelse(level=="gene", "gene_id", "transcript_id"),
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
  counts
}

RPFsCounts <- function(files, txdb, level, bestpsite, readsLen){
  yieldSize <- 10000000
  cnts <- lapply(files, function(f){
    bamfile <- BamFile(file = f, yieldSize = yieldSize)
    pc <- getPsiteCoordinates(bamfile, bestpsite=bestpsite)
    pc.sub <- pc[pc$qwidth %in% readsLen]
    pc.sub <- shiftReadsByFrame(pc.sub, txdb)
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
  features <- switch(level,
                     gene={
                       f <- rep(features, lengths(features$gene_id))
                       mcols(f) <- DataFrame(feature_id=unlist(features$gene_id))
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
  Chr <- sapply(f, foo, n="seqnames")
  Start <- sapply(f, foo, n="start")
  End <- sapply(f, foo, "end")
  Strand <- sapply(f, foo, "strand")
  Length <- sapply(f, function(.ele) sum(.ele$width))
  data.frame(GeneID=GeneID, Chr=Chr, Start=Start, End=End,
             Strand=Strand, Length=Length, stringsAsFactors = FALSE)
}
