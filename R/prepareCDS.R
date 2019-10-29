#' prepare CDS
#' @description prepare CDS library from a TxDb object.
#' @param txdb an TxDb object.
#' @param withUTR including UTR information or not.
#' @return an GRanges object with metadata which include:
#' tx_id: transcript id;
#' tx_name: transcript name;
#' gene_id: gene id;
#' isFirstExonInCDS: is first exon in CDS or not;
#' idFirstExonInCDS: the id for the first exon;
#' isLastExonInCDS: is last exon in CDS or not;
#' wid.cumsu: cumulative sums of number of bases in CDS;
#' internalPos: offset position from 1 base;
#' @importClassesFrom GenomicFeatures TxDb
#' @importFrom GenomicFeatures cdsBy fiveUTRsByTranscript
#' @importFrom AnnotationDbi select
#' @importFrom methods as is
#' @import GenomicRanges
#' @export
#' @examples
#' library(AnnotationDbi)
#' txdb <- loadDb(system.file("extdata", "danRer10.chr1.txdb",
#'                package="ribosomeProfilingQC"))
#' CDS <- prepareCDS(txdb)
#'
prepareCDS <- function(txdb, withUTR = FALSE){
  stopifnot(is(txdb, "TxDb"))
  if(withUTR){
    cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
    utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE)
    utr3 <- threeUTRsByTranscript(txdb, use.names = TRUE)
    CDS <- unlist(cds)
    CDS$cds_id <- NULL
    CDS$cds_name <- NULL
    CDS$feature <- "CDS"
    CDS$tx_name <- rep(names(cds), lengths(cds))
    utrs <- list(UTR5=utr5, UTR3=utr3)
    utrs <- mapply(utrs, names(utrs), FUN=function(.ele, .name){
      UTR <- unlist(.ele)
      UTR$exon_id <- NULL
      UTR$exon_name <- NULL
      UTR$feature <- .name
      UTR$tx_name <- rep(names(.ele), lengths(.ele))
      UTR
    })
    utrs <- unlist(GRangesList(utrs), use.names = FALSE)

    exons <- c(utrs, CDS)
    exons <- sort(exons)
    exons$tx_id <- as.numeric(factor(exons$tx_name, levels = unique(exons$tx_name)))
    exons <- exons[order(exons$tx_id, start(exons)*ifelse(strand(exons)=="-", -1, 1))]

    suppressMessages(id_map <- select(txdb, keys=unique(exons$tx_name),
                                      columns = c("TXNAME", "GENEID", "TXTYPE"), keytype="TXNAME"))
    exons$gene_id <- id_map[match(exons$tx_name, id_map$TXNAME), "GENEID"]
    exons$tx_type <- id_map[match(exons$tx_name, id_map$TXNAME), "TXTYPE"]
    exons$oid <- paste(exons$tx_name, exons$feature)
    isCDS <- exons$feature == "CDS"
    exons$isFirstExonInCDS <- isCDS & !duplicated(exons$oid)
    exons$idFirstExonInCDS <- rep(which(exons$isFirstExonInCDS), table(exons$tx_id)) ## make sure that table is sorted
    exons$isLastExonInCDS <- isCDS & rev(!duplicated(rev(exons$oid)))
    exons$isFirstExonInTx <- !duplicated(exons$tx_id)
    exons$idFirstExonInTx <- rep(which(exons$isFirstExonInTx), table(exons$tx_id))
    exons$isLastExonInTx <- rev(!duplicated(rev(exons$tx_id)))
    exons.wid <- width(exons)
    isUTR5 <- exons$feature == "UTR5"
    CDS.wid.cumsum <- cumsum(exons.wid[!isUTR5])
    exons$wid.cumsum <- 0
    exons.sub <- exons[!isUTR5]
    idFirstExonInCDS <- rep(which(exons.sub$isFirstExonInCDS), table(exons.sub$tx_id))
    exons$wid.cumsum[!isUTR5] <- CDS.wid.cumsum - c(0, CDS.wid.cumsum)[idFirstExonInCDS]
    exons$internalPos <- 0
    exons$internalPos[!(isUTR5 | exons$isFirstExonInCDS)] <-
      exons$wid.cumsum[!(isUTR5 | exons$isLastExonInTx)]
    ## for UTR5
    exons.sub <- rev(exons[isUTR5])
    exons.sub$tx_id <- as.numeric(factor(exons.sub$tx_name, levels = unique(exons.sub$tx_name)))
    exons.sub$isFirstUTR5 <- !duplicated(exons.sub$oid)
    exons.sub$idFirstUTR5 <- rep(which(exons.sub$isFirstUTR5), table(exons.sub$tx_id))
    exons.sub$isLastUTR5 <- rev(!duplicated(rev(exons.sub$oid)))
    exons.sub.wid <- width(exons.sub)
    exons.sub.cumsum <- cumsum(exons.sub.wid)
    exons.sub$wid.cumsum <- exons.sub.cumsum - c(0, exons.sub.cumsum)[exons.sub$idFirstUTR5]
    exons.sub$internalPos <- 0
    exons.sub$internalPos[!exons.sub$isFirstUTR5] <- exons.sub$wid.cumsum[!exons.sub$isLastUTR5]
    exons$wid.cumsum[isUTR5] <-  -1 * rev(exons.sub$wid.cumsum)
    exons$internalPos[isUTR5] <- -1 * rev(exons.sub$internalPos)
    exons$oid <- NULL
    return(exons)
  }
  cds <- cdsBy(txdb, by="tx", use.names = TRUE)
  CDS <- unlist(cds)
  CDS$cds_name <- NULL
  CDS$tx_id <- rep(seq_along(cds), lengths(cds))
  CDS$tx_name <- rep(names(cds), lengths(cds))
  suppressMessages(id_map <- select(txdb, keys=unique(CDS$tx_name),
                                    columns = c("TXNAME", "GENEID", "TXTYPE"), keytype="TXNAME"))
  CDS$gene_id <- id_map[match(CDS$tx_name, id_map$TXNAME), "GENEID"]
  CDS$tx_type <- id_map[match(CDS$tx_name, id_map$TXNAME), "TXTYPE"]
  CDS$isFirstExonInCDS <- !duplicated(CDS$tx_id)
  CDS$idFirstExonInCDS <- rep(which(!duplicated(CDS$tx_id)), lengths(cds))
  CDS$isLastExonInCDS <- rev(!duplicated(rev(CDS$tx_id)))
  CDS.wid <- width(CDS)
  CDS.wid.cumsum <- cumsum(CDS.wid)
  CDS$wid.cumsum <- CDS.wid.cumsum - c(0, CDS.wid.cumsum)[CDS$idFirstExonInCDS]
  CDS$internalPos <- 0
  CDS$internalPos[!CDS$isFirstExonInCDS] <- CDS$wid.cumsum[!CDS$isLastExonInCDS]
  CDS
}
