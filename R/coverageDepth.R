#' extract coverage depth for gene level or transcript level
#' @description Calculate the coverage depth for gene level or transcript level.
#' Coverage for RPFs will be the best P site coverage.
#' Coverage for RNAs will be the coverage for 5'end of reads.
#' @param RPFs Bam file names of RPFs.
#' @param RNAs Bam file names of RNAseq.
#' @param gtf GTF file name for annotation or a TxDb object.
#' @param level transcript or gene level.
#' @param bestpsite P site postion.
#' @param readsLen reads length to keep.
#' @param region annotation region. It could be "cds", "utr5", "utr3",
#' "exon", "transcripts", "feature with extension".
#' @param ext extesion region for "feature with extension".
#' @param ... parameters pass to
#' \link[GenomicFeatures:makeTxDbFromGFF]{makeTxDbFromGFF}
#' @return a list with reads counts.
#' @importFrom methods as is
#' @importFrom GenomicFeatures makeTxDbFromGFF coverageByTranscript transcripts
#' @importFrom Rsamtools BamFile
#' @importFrom S4Vectors DataFrame Rle
#' @importFrom GenomeInfoDb seqnames seqlevels seqlevels<-
#' @importFrom IRanges disjoin
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cvgs <- coverageDepth(RPFs[1], gtf=gtf, level="gene")

coverageDepth <- function(RPFs, RNAs, gtf, level=c("tx", "gene"),
                          bestpsite=13, readsLen=c(28,29), region="cds",
                          ext=5000,
                          ...){
  stopifnot(is.character(gtf))
  level <- match.arg(level)
  region <- region[1]
  if(!region %in% c("cds", "utr5", "utr3", "exon",
                    "transcripts", "feature with extension")){
    stop("region must be cds, utr5, utr3, exon,
         transcripts, feature with extension")
  }
  gtf <- gtf[1]
  cvgs <- list()
  txdb <- NULL
  if(is.character(gtf)){
    suppressWarnings(suppressMessages(txdb <- makeTxDbFromGFF(gtf, ...)))
  }else{
    if(is(gtf, "TxDb")){
      txdb <- gtf
    }
  }
  if(!is(txdb, "TxDb")){
    stop("Can not determine annotations from gtf parameter.")
  }
  if(!missing(RPFs)){
    stopifnot(is.character(RPFs))
    stopifnot(is.numeric(readsLen))
    stopifnot(is.numeric(bestpsite))
    cd <- getCvgs(files = RPFs,
                  txdb = txdb, level = level,
                  bestpsite = bestpsite,
                  readsLen = readsLen,
                  region = region,
                  ext = ext)
    cvgs[["RPFs"]] <- cd
  }
  if(!missing(RNAs)){
    stopifnot(is.character(RNAs))
    cd <- getCvgs(files = RNAs,
                  txdb = txdb, level = level,
                  region = region,
                  ext = ext)
    cvgs[["mRNA"]] <- cd
  }
  cvgs
}

getCvgs <- function(files, txdb, level, bestpsite, readsLen, region, ext){
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
  names(reads) <- basename(sub(".bam", "", files, ignore.case = TRUE))
  features <- switch(region,
    cds=cds(txdb, columns=c("gene_id", "tx_name")),
    utr5={
      f <- fiveUTRsByTranscript(txdb, use.name=TRUE)
      fs <- unlist(f)
      mcols(fs) <- DataFrame(tx_name=rep(names(f), lengths(f)))
      suppressMessages(id_map <-
                         select(txdb, keys=unique(fs$tx_name),
                                columns = c("TXNAME", "GENEID"),
                                keytype="TXNAME"))
      fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME), "GENEID"]
      fs
    },
    utr3={
      f <- threeUTRsByTranscript(txdb, use.name=TRUE)
      fs <- unlist(f)
      mcols(fs) <- DataFrame(tx_name=rep(names(f), lengths(f)))
      suppressMessages(id_map <-
                         select(txdb, keys=unique(fs$tx_name),
                                columns = c("TXNAME", "GENEID"),
                                keytype="TXNAME"))
      fs$gene_id <- id_map[match(fs$tx_name, id_map$TXNAME), "GENEID"]
      fs
    },
    exon=exons(txdb, columns=c("gene_id", "tx_name")),
    transcripts=transcripts(txdb, columns=c("gene_id", "tx_name")),
    "feature with extension"={
      f <- transcripts(txdb, columns=c("gene_id", "tx_name"))
      e <- exons(txdb, columns=c("gene_id", "tx_name"))
      tx_name <- e$tx_name
      e <- rep(e, lengths(tx_name))
      e$tx_name <- unlist(tx_name)
      CDS <- cds(txdb, columns=c("gene_id", "tx_name"))
      tx_name <- CDS$tx_name
      CDS <- rep(CDS, lengths(tx_name))
      CDS$tx_name <- unlist(tx_name)
      ## avoid overlaps of extension region with feature region
      getExt <- function(ext, f, e, start=TRUE){
        exts <- flank(f, ext, start=start)
        diffs <- setdiff(exts, f)
        ols <- findOverlaps(diffs, exts, type = "within")
        es <- diffs[queryHits(ols)]
        mcols(es) <- mcols(exts[subjectHits(ols)])
        dists <- distance(es, f[match(es$tx_name, f$tx_name)])
        es <- es[dists==0]
        ## merge to exon
        ol <- findOverlaps(es, e, maxgap = 1)
        ol <- ol[es[queryHits(ol)]$tx_name == e[subjectHits(ol)]$tx_name]
        es <- es[queryHits(ol)]
        e1 <- e[subjectHits(ol)]
        dist <- distance(es, e1)
        es <- es[dist==0]
        e1 <- e1[dist==0]
        ol <- ol[dist==0]
        start(es) <- ifelse(start(es)<start(e1), start(es), start(e1))
        end(es) <- ifelse(end(es)<end(e1), end(e1), end(es))
        e[subjectHits(ol)] <- es
        e
      }
      exStart <- getExt(ext, f, e, start=TRUE)
      exEnd <- getExt(ext, f, e, start=FALSE)
      f_ex <- c(f, CDS, exStart, exEnd)
      f_ex <- split(f_ex, f_ex$tx_name)
      f_ex <- disjoin(f_ex)
      f_ex <- unlist(f_ex)
      mcols(f_ex) <- mcols(f)[match(names(f_ex), f$tx_name), ]
      f_ex
    }
  )
  features <- switch(level,
                     gene={
                       f <- rep(features, lengths(features$gene_id))
                       mcols(f) <-
                         DataFrame(feature_id=unlist(features$gene_id))
                       f[!is.na(f$feature_id)]
                     },
                     tx={
                       f <- rep(features, lengths(features$tx_name))
                       mcols(f) <-
                         DataFrame(feature_id=unlist(features$tx_name))
                       f[!is.na(f$feature_id)]
                     })
  cvgs <- lapply(reads, coverage)
  seqs <- lapply(cvgs, function(.ele)
    names(.ele)[vapply(.ele, sum, FUN.VALUE = 0)>0])
  seqs <- unique(unlist(seqs))
  seqs <- intersect(seqlevels(features), seqs)
  features <- features[as.character(seqnames(features)) %in% seqs]
  seqlevels(features) <- seqlevels(features)[seqlevels(features) %in% seqs]
  features1 <- disjoin(features, with.revmap=TRUE)
  f <- rep(features1, lengths(features1$revmap))
  f$feature_id <-
    unlist(lapply(features1$revmap, function(.ele) features$feature_id[.ele]))
  f$revmap <- NULL
  f <- f[!duplicated(f) | !duplicated(f$feature_id)]
  ## make sure that all features are sorted by start position.
  f <- f[order(f$feature_id, as.character(seqnames(f)),
               ifelse(as.character(strand(f))=="-", -1, 1)*start(f))]
  f <- split(f, f$feature_id)
  coverages <- lapply(reads, coverageByTranscript, transcripts=f,
                      ignore.strand=ignore.strand)
  cd <- cvgd(coverage=coverages, granges=f)
  return(cd)
}
