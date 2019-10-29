#' extract coverage depth for gene level or transcript level
#' @description Calculate the coverage depth for gene level or transcript level.
#' Coverage for RPFs will be the best P site coverage.
#' Coverage for RNAs will be the coverage for 5'end of reads.
#' @param RPFs Bam file names of RPFs.
#' @param RNAs Bam file names of RNAseq.
#' @param gtf GTF file name for annotation.
#' @param level transcript or gene level.
#' @param bestpsite P site postion.
#' @param readsLen reads length to keep.
#' @param region annotation region. It could be "cds", "utr5", "utr3" or "exon".
#' @param ... parameters pass to \link[Rsubread:featureCounts]{featureCounts}
#' @return a list with reads counts.
#' @importFrom methods as is
#' @importFrom Rsubread featureCounts
#' @importFrom GenomicFeatures makeTxDbFromGFF coverageByTranscript
#' @importFrom Rsamtools BamFile
#' @importFrom S4Vectors DataFrame Rle
#' @importFrom GenomeInfoDb seqnames seqlevels seqlevels<-
#' @importFrom IRanges disjoin
#' @export
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cvgs <- coverageDepth(RPFs, RNAs, gtf, level="tx")
#' }

coverageDepth <- function(RPFs, RNAs, gtf, level=c("tx", "gene"),
                          bestpsite=13, readsLen=c(28,29), region="cds",
                          ...){
  stopifnot(is.character(gtf))
  level <- match.arg(level)
  region <- region[1]
  if(!region %in% c("cds", "utr5", "utr3", "exon")){
    stop("region must be cds, utr5, utr3, or exon")
  }
  gtf <- gtf[1]
  cvgs <- list()
  suppressWarnings(suppressMessages(txdb <- makeTxDbFromGFF(gtf)))
  if(!missing(RPFs)){
    stopifnot(is.character(RPFs))
    stopifnot(is.numeric(readsLen))
    stopifnot(is.numeric(bestpsite))
    cvgs[["RPFs"]] <- getCvgs(files = RPFs,
                               txdb = txdb, level = level,
                               bestpsite = bestpsite,
                               readsLen = readsLen,
                              region = region)
  }
  if(!missing(RNAs)){
    stopifnot(is.character(RNAs))
    cvgs[["mRNA"]] <- getCvgs(files = RNAs,
                              txdb = txdb, level = level,
                              region = region)
  }
  cvgs
}

getCvgs <- function(files, txdb, level, bestpsite, readsLen, region){
  yieldSize <- 10000000
  if(!missing(bestpsite)){##RPFs
    ignore.strand <- FALSE
    reads <- lapply(files, function(f){
      bamfile <- BamFile(file = f, yieldSize = yieldSize)
      pc <- getPsiteCoordinates(bamfile, bestpsite=bestpsite)
      pc.sub <- pc[pc$qwidth %in% readsLen]
      pc.sub <- shiftReadsByFrame(pc.sub, txdb)
    })
  }else{##mRNA
    ignore.strand <- TRUE
    reads <- lapply(files, function(f){
      bamfile <- BamFile(file = f, yieldSize = yieldSize)
      pc <- getPsiteCoordinates(bamfile, bestpsite=1)
    })
  }
  names(reads) <- basename(sub(".bam", "", files))
  features <- switch(region,
    cds=cds(txdb, columns=c("gene_id", "tx_name")),
    utr5={
      f <- fiveUTRsByTranscript(txdb, use.name=TRUE)
      fs <- unlist(f)
      mcols(fs) <- DataFrame(tx_name=rep(names(f), lengths(f)))
      suppressMessages(id_map <- select(txdb, keys=unique(fs$tx_name),
                                        columns = c("TXNAME", "GENEID"), keytype="TXNAME"))
      fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME), "GENEID"]
      fs
    },
    utr3={
      f <- threeUTRsByTranscript(txdb, use.name=TRUE)
      fs <- unlist(f)
      mcols(fs) <- DataFrame(tx_name=rep(names(f), lengths(f)))
      suppressMessages(id_map <- select(txdb, keys=unique(fs$tx_name),
                                        columns = c("TXNAME", "GENEID"), keytype="TXNAME"))
      fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME), "GENEID"]
      fs
    },
    exon=exons(txdb, columns=c("gene_id", "tx_name"))
  )
  features <- switch(level,
                     gene={
                       f <- rep(features, lengths(features$gene_id))
                       mcols(f) <- DataFrame(feature_id=unlist(features$gene_id))
                       f[!is.na(f$feature_id)]
                     },
                     tx={
                       f <- rep(features, lengths(features$tx_name))
                       mcols(f) <- DataFrame(feature_id=unlist(features$tx_name))
                       f[!is.na(f$feature_id)]
                     })
  cvgs <- lapply(reads, coverage)
  seqs <- lapply(cvgs, function(.ele) names(.ele)[sapply(.ele, sum)>0])
  seqs <- unique(unlist(seqs))
  seqs <- intersect(seqlevels(features), seqs)
  features <- features[as.character(seqnames(features)) %in% seqs]
  seqlevels(features) <- seqlevels(features)[seqlevels(features) %in% seqs]
  features1 <- disjoin(features, with.revmap=TRUE)
  f <- rep(features1, lengths(features1$revmap))
  f$feature_id <- unlist(lapply(features1$revmap, function(.ele) features$feature_id[.ele]))
  f$revmap <- NULL
  f <- f[!duplicated(f) | !duplicated(f$feature_id)]
  f <- f[order(f$feature_id, as.character(seqnames(f)),
               ifelse(as.character(strand(f))=="-", -1, 1)*start(f))] ## make sure that all features are sorted by start position.
  f <- split(f, f$feature_id)
  return(lapply(reads, coverageByTranscript, transcripts=f, ignore.strand=ignore.strand))
}
