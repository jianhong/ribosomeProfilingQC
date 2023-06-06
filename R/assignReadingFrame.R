#' Assign reading frame
#' @description Set reading frame for each reads in CDS region to frame0,
#' frame1 and frame2.
#' @param reads Output of \link{getPsiteCoordinates}
#' @param CDS Output of \link{prepareCDS}
#' @param txdb A TxDb object. If it is set, assign reading frame for all reads.
#' Default missing, only assign rading frame for reads in CDS.
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @return An GRanges object of reads with reading frame information.
#' @import GenomicRanges
#' @importFrom methods as is
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @export
#' @examples
#' library(Rsamtools)
#' bamfilename <- system.file("extdata", "RPF.WT.1.bam",
#'                            package="ribosomeProfilingQC")
#' yieldSize <- 10000000
#' bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
#' pc <- getPsiteCoordinates(bamfile, bestpsite=13)
#' pc.sub <- pc[pc$qwidth %in% c(29, 30)]
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
#' pc.sub <- assignReadingFrame(pc.sub, CDS)

assignReadingFrame <- function(reads, CDS, txdb, ignore.seqlevelsStyle=FALSE){
  if(!missing(txdb)){
    stopifnot(is(txdb, "TxDb"))
    CDS <- prepareCDS(txdb, withUTR=TRUE)
  }
  stopifnot(is(CDS, "GRanges"))
  if(length(CDS$internalPos)!=length(CDS) ||
     length(CDS$isLastExonInCDS)!=length(CDS) ||
     length(CDS$tx_name)!=length(CDS) ||
     length(CDS$gene_id)!=length(CDS)){
    stop("CDS must be output of prepareCDS")
  }
  stopifnot(is(reads, "GRanges"))
  reads <-fixSeqlevelsStyle(reads, CDS, ignore.seqlevelsStyle)
  reads$readingFrame <- NA
  reads$position <- NA
  reads$offsetPercentage <- NA
  reads$posToStop <- NA
  reads$tx_name <- NA
  reads$gene_id <- NA
  CDS.last <- CDS[CDS$isLastExonInCDS]
  names(CDS.last) <- CDS.last$tx_name
  ## check reading frame in CDS
  ol <- findOverlaps(reads, CDS, ignore.strand=FALSE)
  sh <- CDS[subjectHits(ol)]
  qh <- reads[queryHits(ol)]
  dist <- distance(qh, promoters(sh, upstream = 1, downstream = 0),
                   ignore.strand=FALSE)
  if(!missing(txdb)){
    sh.rev <- switch.strand(sh)
    p <- promoters(sh.rev, upstream = 1, downstream = 0)
    p <- switch.strand(p)
    dist2 <- -1 * distance(qh, p, ignore.strand=FALSE)
    dist <- ifelse(sh$wid.cumsum<0, dist2, dist)
  }
  pos <- dist + sh$internalPos
  cvg <- pos/CDS.last[sh$tx_name]$wid.cumsum
  toStop <- CDS.last[sh$tx_name]$wid.cumsum - pos - 3 # stop codon
  tx_name <- sh$tx_name
  gene_id <- sh$gene_id
  pos <- data.frame(id=queryHits(ol), frame=pos %% 3, pos=pos, cvg=cvg,
                    toStop=toStop, tx_name=tx_name, gene_id=gene_id,
                    stringsAsFactors=FALSE)
  pos <- pos[order(pos$id, pos$frame), ]
  pos <- pos[!duplicated(pos$id), ]
  reads$readingFrame[pos$id] <- pos$frame
  reads$position[pos$id] <- pos$pos
  reads$offsetPercentage[pos$id] <- pos$cvg
  reads$posToStop[pos$id] <- pos$toStop
  reads$tx_name[pos$id] <- pos$tx_name
  reads$gene_id[pos$id] <- pos$gene_id
  reads
}
