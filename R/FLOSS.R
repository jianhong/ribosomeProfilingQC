#' Fragment Length Organization Similarity Score (FLOSS)
#' @description The FLOSS will be calculated from a histogram of
#' read lengths for footprints on a transcript or reading frame.
#' @references 1: Ingolia NT, Brar GA, Stern-Ginossar N, Harris MS,
#' Talhouarne GJ, Jackson SE, Wills MR, Weissman JS. Ribosome profiling
#' reveals pervasive translation outside
#' of annotated protein-coding genes. Cell Rep. 2014 Sep 11;8(5):1365-79. doi:
#' 10.1016/j.celrep.2014.07.045. Epub 2014 Aug 21. PubMed PMID: 25159147; PubMed
#' Central PMCID: PMC4216110.
#' @param reads Output of \link{getPsiteCoordinates}
#' @param ref Refercence id list. If level is set to tx, the id should be
#' transcript names. If level is set to gene, the id should be gene id.
#' @param CDS Output of \link{prepareCDS}
#' @param readLengths Read length used for calculation
#' @param level Transcript or gene level
#' @param draw Plot FLOSS vs total reads or not.
#' @return A data frame with colnames as id, FLOSS, totalReads,
#' wilcox.test.pval, cook's distance.
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stats cooks.distance lm wilcox.test
#' @importFrom ggplot2 geom_point geom_smooth scale_x_log10 scale_y_log10
#' annotation_logticks
#' @importFrom ggrepel geom_label_repel
#' @importFrom scales trans_breaks trans_format math_format
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.WT.1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' pc <- getPsiteCoordinates(bamfile, bestpsite=13)
#' #library(GenomicFeatures)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' #txdb <- makeTxDbFromGFF(system.file("extdata",
#'  #         "Danio_rerio.GRCz10.91.chr1.gtf.gz",
#'  #         package="ribosomeProfilingQC"),
#'  #         organism = "Danio rerio",
#'  #         chrominfo = seqinfo(Drerio)["chr1"],
#'  #         taxonomyId = 7955)
#' #CDS <- prepareCDS(txdb)
#' CDS <- readRDS(system.file("extdata", "CDS.rds",
#'                            package="ribosomeProfilingQC"))
#' set.seed(123)
#' ref <- sample(unique(CDS$gene_id), 100)
#' fl <- FLOSS(pc, ref, CDS, level="gene")

FLOSS <- function(reads, ref, CDS, readLengths=c(26:34),
                  level=c("tx", "gene"), draw=FALSE){
  stopifnot(is(reads, "GRanges"))
  stopifnot(length(reads$qwidth)==length(reads))
  level <- match.arg(level)
  level <- c("tx"="tx_name", "gene"="gene_id")[level]
  stopifnot(is(CDS, "GRanges"))
  if(length(CDS$tx_name)!=length(CDS) ||
     length(CDS$gene_id)!=length(CDS)){
    stop("CDS must be output of prepareCDS")
  }
  if(length(intersect(seqlevelsStyle(reads), seqlevels(CDS)))==0){
    seqlevelsStyle(reads) <- seqlevelsStyle(CDS)[1]
  }
  stopifnot(is.numeric(readLengths))
  readLengths <- as.integer(readLengths)
  reads <- reads[reads$qwidth>=min(readLengths) &
                   reads$qwidth<=max(readLengths)]
  if(length(reads)<1){
    stop("No reads is in given readLengths.")
  }
  ol <- findOverlaps(CDS, reads)
  qh <- mcols(CDS[queryHits(ol)])[, level]
  sh <- reads[subjectHits(ol)]
  sh$ID <- qh
  sh <- sh[(!duplicated(sh)) | (!duplicated(qh))]
  cnt <- data.frame(len=sh$qwidth, ID=sh$ID)
  cnt.ref <- cnt[cnt$ID %in% ref, "len"]
  cnt.tab <- table(cnt)
  cnt.ref <- table(cnt.ref)
  cnt.ref <- cnt.ref[match(as.character(readLengths), names(cnt.ref))]
  cnt.ref[is.na(cnt.ref)] <- 0
  names(cnt.ref) <- as.character(readLengths)
  cnt.ref.sum <- sum(cnt.ref)
  if(cnt.ref.sum<1){
    stop("No reads is in ref.")
  }
  f.ref <- cnt.ref/cnt.ref.sum
  cnt.tab <- cnt.tab[match(as.character(readLengths), rownames(cnt.tab)), ]
  rownames(cnt.tab) <- as.character(readLengths)
  cnt.tab[is.na(cnt.tab)] <- 0
  cnt.tab <- cnt.tab[, colSums(cnt.tab)>0, drop=FALSE]
  cnt.tab.sum <- colSums(cnt.tab)
  f.tab <- t(cnt.tab)/cnt.tab.sum

  fl <- apply(f.tab, 1, function(.ele) .ele-f.ref)
  fl <- .5 * colSums(abs(fl))
  pval <- apply(f.tab, 1,
                function(.ele){
                  .id <- (.ele - f.ref)!=0
                  if(sum(.id)==0){
                    1
                  }else{
                    wilcox.test(x = .ele[.id],
                                y = f.ref[.id], paired = TRUE)$p.value
                  }
                })
  id <- names(fl)
  df <- data.frame(id=id, FLOSS=fl, totalReads=cnt.tab.sum[id],
                   p.value=pval[id], stringsAsFactors = FALSE)
  fit.data <- log2(df[df$FLOSS!=0, c("totalReads", "FLOSS")])
  fit <- lm(FLOSS~0+totalReads,
            log2(df[df$FLOSS!=0, c("totalReads", "FLOSS")]))
  cooksd <- cooks.distance(fit)
  df$cooks.distance <- 0
  df$cooks.distance[df$FLOSS!=0] <- cooksd
  if(draw){
    ggdf <- df[df$FLOSS!=0,
               c("id", "totalReads", "FLOSS", "cooks.distance"),
               drop=FALSE]
    ggdf$highlight <- ggdf$cooks.distance>4*mean(cooksd, na.rm=TRUE)
    colnames(ggdf) <- c("ggid", "ggx", "ggy", "ggcd", "gghighlight")
    ggid <- ggx <- ggy <- ggcd <- gghighlight <- .x <- NULL
    plot <- ggplot(ggdf, aes(x=ggx, y=ggy)) +
      geom_point(pch = 16, cex = .5) +
      geom_smooth(method = "lm", formula = y~0+x, color="red", lty=3) +
      xlab("total reads") +
      ylab("FLOSS") +
      geom_label_repel(data = subset(ggdf, gghighlight), aes(label=ggid)) +
      geom_point(data = subset(ggdf, gghighlight), pch=21, color="red") +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      annotation_logticks() +
      theme_classic()
    print(plot)
  }
  df <- cbind(df,
              matrix(as.numeric(f.tab[id, ]),
                     ncol=ncol(f.tab),
                     dimnames = list(ID=id, len=colnames(f.tab))))
  df.ref=data.frame(id="ref", FLOSS=0, totalReads=cnt.ref.sum,
                 p.value=1, cooks.distance=0,
                 t(unlist(as.list(f.ref))), row.names = "ref",
                 stringsAsFactors = FALSE, check.names = FALSE)
  df <- rbind(df.ref, df)
  return(df)
}
