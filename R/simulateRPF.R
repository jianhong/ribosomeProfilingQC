#' simulate function
#' @description simulate the RPFs reads in CDS, 5'UTR and 3'UTR
#' @param txdb TxDb object
#' @param outPath output folder for the bam files
#' @param genome BSgenome object
#' @param samples Total samples to simulate.
#' @param group1,group2 numeric to index the sample groups.
#' @param readsPerSample Total reads number per sample.
#' @param readsLen Reads length, default 100bp.
#' @param psite P-site position. default 13.
#' @param frame0,frame1,frame2 Percentage of reads distribution in
#' frame0, frame1 and frame2
#' @param DEregions The regions with differential reads in exon,
#' utr5 and utr3.
#' @param size dispersion parameter. Must be strictly positive, need not be integer.
#' @param sd standard deviations.
#' @param includeReadsSeq logical. Include reads sequence or not.
#' @return an invisible list of GAlignments.
#' @importFrom rtracklayer export
#' @importFrom GenomicAlignments sequenceLayer cigar GAlignments
#' @importFrom BSgenome getSeq
#' @importFrom S4Vectors runLength
#' @importFrom stats rnbinom rnorm runif
#' @importFrom Biostrings DNAStringSet
#' @importFrom GenomeInfoDb seqlengths `seqlengths<-`
#' @export
#' @examples
#' library(GenomicFeatures)
#' txdb <- loadDb(system.file("extdata", "danRer10.chr1.txdb",
#'                package="ribosomeProfilingQC"))
#' simulateRPF(txdb, samples=1, readsPerSample = 1e3)
#' \dontrun{
#' cds <- prepareCDS(txdb, withUTR = TRUE)
#' cds <- cds[width(cds)>100]
#' DEregions <- cds[sample(seq_along(cds), 500)]
#' simulateRPF(txdb, samples=6, readsPerSample = 1e5, DEregions=DEregions)
#' }

simulateRPF <- function(txdb, outPath, genome, samples = 6,
                        group1 = 1:3, group2 = 4:6,
                        readsPerSample = 1e6, readsLen = 28,
                        psite = 13,
                        frame0=.90, frame1=.05, frame2=.05,
                        DEregions=GRanges(),
                        size = 2, sd = .05,
                        includeReadsSeq = FALSE){
  cds <- prepareCDS(txdb, withUTR = TRUE)
  stopifnot(sum(c(frame0, frame1, frame2))==1)
  stopifnot(psite<readsLen)
  stopifnot(size>0)
  mis_genome <- missing(genome)
  if(!mis_genome) stopifnot(is(genome, "BSgenome"))
  stopifnot(is(DEregions, "GRanges"))
  ol_DE <- findOverlaps(DEregions, cds, type = "within")
  if(any(!seq_along(DEregions) %in% queryHits(ol_DE))){
    stop('Not all DEregions within exon regions.')
  }
  if(is.numeric(samples)){
    samples <- paste0("sample", seq.int(samples[1]))
  }
  stopifnot(is(group1, "numeric"))
  stopifnot(is(group2, "numeric"))
  if(length(DEregions)>0){
    stopifnot(length(group1)==length(group2))
    if(any(!group1 %in% seq_along(samples)) ||
       any(!group2 %in% seq_along(samples))){
      stop("index of group1 and group2 must within sample numbers.")
    }
    if(!all(seq_along(samples) %in% c(group1, group2))){
      stop("index of group1 and group2 must cover all samples.")
    }
  }
  ## set reads in CDS
  CDS <- cds[cds$feature=="CDS"]
  ## get all cds cumsum
  CDS.length <- CDS$wid.cumsum[CDS$isLastExonInCDS]
  names(CDS.length) <- CDS$tx_name[CDS$isLastExonInCDS]
  ## check all cds length mod by 3
  check_CDS <- (CDS.length %% 3) != 0
  if(any(check_CDS)){
    warning("There are ", sum(check_CDS), "(",
            round(100*(sum(check_CDS)/length(check_CDS)), 2), "%) ",
            "not mudulo3 CDSs (which length can not by divided by 3):",
            paste(names(CDS.length)[check_CDS], collapse = ", "))
  }
  CDS.length.3 <- CDS.length %/% 3
  ### utr5 and utr3
  u5 <- cds[cds$isFirstExonInTx & !cds$isFirstExonInCDS]
  u5Len <- rep(0, length(unique(cds$tx_name)))
  names(u5Len) <- unique(cds$tx_name)
  u5Len[u5$tx_name] <- abs(u5$wid.cumsum)
  fullLen <- split(cds, cds$tx_name)
  fullLen <- sapply(fullLen, function(.ele) sum(width(.ele)))
  fullLen <- fullLen - readsLen

  ## random model, negative binomial distribution
  n <- length(CDS.length.3)
  x <- rnbinom(n, mu=readsPerSample/n, size = size)
  ## normal distribution for samples
  y <- lapply(x, function(mu) round(rnorm(n=length(samples), mean = mu, sd=sd*mu)))

  ## assign position
  p <- mapply(CDS.length.3, y, u5Len[names(CDS.length.3)], fullLen[names(CDS.length.3)],
              FUN = function(m, n, ups, dws){
    lapply(n, function(.n){
      .min <- max(0, ceiling((psite-1-ups)/3))## if there is no utr5, it should leave the space for p0 to psite
      .max <- min(m-1, max(1, floor((dws+psite-1-ups)/3))) ## if there is no utr3, it should leave the space from psite to end
      if(.min>.max) return(NULL)
      3*round(runif(n=.n,
                    min=.min,
                    max=.max))
    })
  }, SIMPLIFY = FALSE)


  ## swap list by sample in level 1
  p <- lapply(seq_along(samples), function(i){
    lapply(p, function(.ele) .ele[[i]])
  })
  names(p) <- samples
  ## sample for frame0, 1, 2
  ## and then thift to p0
  p <- lapply(p, function(.sample){
    pos <- GRanges(rep(names(.sample), lengths(.sample)),
                       IRanges(unlist(.sample), width = 1))
    offset <- sample(x = c(0, 1, 2),
                     size = length(pos),
                     replace = TRUE,
                     prob = c(frame0, frame1, frame2))
    pos <- shift(pos, shift = offset-psite+1)
  })

  if(length(samples)>1){
    cds_de <- cds[subjectHits(ol_DE)]
    cds_de <- split(cds_de, cds_de$feature)
    ## add de 5'UTR and 3'UTRs
    ## sample the reads from CDS region and shift it into 5'UTR or 3'UTR for up group
    ## all psite is in CDS now
    ## 1. get the length of 5'UTR and 3'UTR > readsLen
    ## 2. split it into two group up and down
    ## 3. get total reads for each txs, and reditribute the reads in 5'UTR (or 3'UTR) + CDS
    if(length(cds_de$UTR3)>0 || length(cds_de$UTR5)>0){
      de_evt <- c(cds_de$UTR3, cds_de$UTR5)
      execpt <- width(de_evt)< readsLen
      if(any(execpt)){
        message("There are ", sum(execpt), " DE UTRs width smaller than readsLen.")
      }
      de_evt <- de_evt[!execpt]
      pp <- unlist(GRangesList(p))
      pp$group <- rep(names(p), lengths(p))
      tx_avg_exp <- lengths(split(pp, seqnames(pp)))
      tx_avg_full <- tx_avg_exp/fullLen[names(tx_avg_exp)]
      de_times <- sign(rnorm(length(de_evt), sd = sd)) *
        (rnbinom(length(de_evt), size=size, prob = .5) + log2(1.5)) +
        rnorm(length(de_evt), sd = sd)
      ## low value
      de_aliq <- abs(de_times) + 1
      de_aliq <- tx_avg_full[de_evt$tx_name]*width(de_evt)/de_aliq
      ## high value
      de_high <- de_aliq * abs(de_times)
      de_to_move <- round(de_aliq + de_high)
      de_aliq <- round(de_aliq)
      de_high <- de_to_move - de_aliq
      ## sample which one to move
      pp <- pp[order(seqnames(pp))]
      pp$id <- paste(as.character(seqnames(pp)),
                     unlist(lapply(runLength(seqnames(pp)), seq.int)))
      de_to_mv_id <- mapply(tx_avg_exp[de_evt$tx_name], de_to_move,
                            FUN=function(n, size){
                              sample.int(n, size)
                            }, SIMPLIFY = FALSE)
      de_to_mv_id <- paste(rep(de_evt$tx_name, de_to_move),
                           unlist(de_to_mv_id))
      ## reasign group
      de_to_mv_gp <- mapply(de_high, de_aliq, de_to_move, sign(de_times),
                            FUN=function(a, b, c, d){
                              x <- samples[c(group2, group1)]
                              if(d<0) x <- samples[c(group1, group2)]
                              sample(x, size=c, replace=TRUE,
                                     prob=rep(c(a, b), each=length(group1)))
                            }, SIMPLIFY = FALSE)
      de_to_mv_gp <- unlist(de_to_mv_gp)
      pp[match(de_to_mv_id, pp$id)]$group <- de_to_mv_gp
      ## move to where?
      ## 3'UTR: internalPos + width
      ## 5'UTR: internalPos - width
      de_to_mv_pos <- mapply(de_to_move,
                             de_evt$internalPos,
                             width(de_evt),
                             de_evt$isLastExonInTx,
                             de_evt$isFirstExonInTx,
                             de_evt$feature,
                             FUN = function(n, m, w, isLast, isFirst, feature){
                               if(isLast) w <- w - readsLen + psite
                               s <- sample.int(w, size = n, replace = TRUE)
                               if(feature=="UTR3"){
                                 return(m + s)
                               }else{
                                 return(m - s)
                               }
                             }, SIMPLIFY = FALSE)
      de_to_mv_pos <- unlist(de_to_mv_pos)
      ranges(pp[match(de_to_mv_id, pp$id)]) <- IRanges(de_to_mv_pos, width = 1)
      gp <- pp$group
      mcols(pp) <- NULL
      pp <- unname(pp)
      p <- split(pp, gp)
    }
    ## add de exons
    ## move the reads from down group to up group
    if(length(cds_de$CDS)>0){
      de_evt <- cds_de$CDS
      de_evt_rg <- GRanges(de_evt$tx_name, IRanges(de_evt$internalPos+1, de_evt$wid.cumsum))
      pp <- unlist(GRangesList(p))
      pp$group <- rep(names(p), lengths(p))
      pp$gp <- pp$group %in% samples[group1]
      pp$id <- seq_along(pp)
      de_evt_ol <- suppressWarnings(findOverlaps(de_evt_rg, pp))
      execpt <- !seq_along(de_evt_rg) %in% queryHits(de_evt_ol)
      if(any(execpt)){
        message("There are ", sum(execpt), " DE GRanges does not contain any simulated reads.")
      }
      de_evt <- de_evt[!execpt]
      de_evt_rg <- de_evt_rg[!execpt]
      de_cnt <- suppressWarnings(countOverlaps(de_evt_rg, pp)/2)
      de_times <- sign(rnorm(length(de_cnt), sd = sd)) *
        (rnbinom(length(de_cnt), size=size, prob = .5) + log2(1.5)) +
        rnorm(length(de_cnt), sd = sd)
      de_aliq <- abs(de_times) + 1
      de_aliq <- 2*de_cnt/de_aliq
      de_to_move <- round((de_aliq*abs(de_times)-de_cnt)*sign(de_times))
      de_evt_ol <- suppressWarnings(findOverlaps(de_evt_rg, pp))
      de_evt_reads <- split(pp[subjectHits(de_evt_ol)], queryHits(de_evt_ol))
      de_to_move <- de_to_move[as.numeric(names(de_evt_reads))]
      notsig <- abs(de_to_move) > lengths(de_evt_reads) | de_to_move == 0
      if(any(notsig)){
        message("There are ", sum(notsig), " DE GRanges does not contain enough simulated reads to show difference.")
      }
      renameM <- samples[c(group1, group2)]
      names(renameM) <- samples[c(group2, group1)]
      de_evt_reads <- de_evt_reads[!notsig]
      de_evt_reads <- mapply(de_to_move[!notsig], de_evt_reads,
                             FUN=function(.n, .reads){
                               .size <- abs(.n)
                               x <- seq_along(.reads)
                               x <- if(sign(.n)>0){
                                 x[!.reads$gp]
                               }else{
                                 x[.reads$gp]
                               }
                               if(length(x)>0){
                                 .id <- sample(x = x, size = .size, replace = TRUE)
                                 .id <- unique(.id)
                                 .reads$group[.id] <- renameM[.reads$group[.id]]
                               }
                               .reads
                             })
      de_evt_reads <- unlist(GRangesList(de_evt_reads))
      pp[match(de_evt_reads$id, pp$id)]$group <- de_evt_reads$group
      gp <- pp$group
      mcols(pp) <- NULL
      pp <- unname(pp)
      p <- split(pp, gp)
    }
  }


  ## add 5'UTR and 3'UTR info
  p <- lapply(p, function(.sample){
    .sample <- shift(.sample, shift = u5Len[as.character(seqnames(.sample))])
    start(.sample)[start(.sample)<0] <- 0 ## no need, but keep
    .sample <- shift(.sample, shift = 1) ## shift 1, why not move to assign position step?
    .sample[start(.sample)>fullLen[as.character(seqnames(.sample))]] <-
      shift(.sample[start(.sample)>fullLen[as.character(seqnames(.sample))]],
            shift = -3)
    start(.sample[start(.sample)>fullLen[as.character(seqnames(.sample))]]) <-
      fullLen[as.character(seqnames(.sample))[
        start(.sample)>fullLen[as.character(seqnames(.sample))]]] ## no need, but keep
    .sample
  })

  #relative2realLoc
  cds.wid.cumsum <- cumsum(width(cds))
  cds$w_cumsum <- cds.wid.cumsum - c(0, cds.wid.cumsum)[cds$idFirstExonInTx]
  cds$i_cumsum <- 0
  cds$i_cumsum[!cds$isFirstExonInTx] <- cds$w_cumsum[!cds$isLastExonInTx]
  cds$i_id <- seq_along(cds)
  cds_rel <- GRanges(cds$tx_name, IRanges(cds$i_cumsum+1, cds$w_cumsum))

  gl <- lapply(p, function(.sample){
    ## convert p to real genomic locations
    getPosInReal <- function(.re){
      ol <- findOverlaps(.re, cds_rel)
      stopifnot(all(queryHits(ol)==seq_along(.re)))
      .pos <- cds[subjectHits(ol)]
      .shift <- start(.re[queryHits(ol)]) - .pos$i_cumsum - 1
      .pos <- promoters(.pos, upstream = 0, downstream = 1)
      .pos <- shift(.pos,
                    shift = ifelse(as.character(strand(.pos))=="-",
                                   -1*.shift,
                                   .shift))
    }
    pos_start <- getPosInReal(.sample)
    pos_stop <- getPosInReal(shift(.sample, shift = readsLen-1))
    ## prepare cigar
    ## most reads only in exon: M28
    ## reads in multple exons: MxNyMz...
    ## i_id_range is the start and end position or the reads
    sw <- start(pos_start)>start(pos_stop)
    st <- ifelse(sw, start(pos_stop), start(pos_start))
    se <- ifelse(sw, start(pos_start), start(pos_stop))
    i_id_range <- rbind(st, se)
    stopifnot(all(i_id_range[1, ] < i_id_range[2, ]))
    i_mult_id <- pos_start$i_id != pos_stop$i_id & ## different feature
      pos_start$exon_rank != pos_stop$exon_rank ## different exon
    i_id_range <- as.list(as.data.frame(i_id_range[, i_mult_id]))
    ## insert covered intron into i_id_range
    ss <- start(cds)
    es <- end(cds)
    std <- as.character(strand(cds))=="-"
    i_id_range <- mapply(i_id_range,
                         pos_start$i_id[i_mult_id],
                         pos_stop$i_id[i_mult_id],
                         FUN=function(ir, id1, id2){
                           .ir <- c(ss[id1:id2], es[id1:id2])
                           .ir <- .ir[.ir>ir[1] & .ir<ir[2]]
                           if(std[id1]){
                             .ir <- c(.ir, es[id2][es[id2] %in% ir], ss[id1][ss[id1] %in% ir])
                           }else{
                             .ir <- c(.ir, es[id1][es[id1] %in% ir], ss[id2][ss[id2] %in% ir])
                           }
                           sort(c(ir, .ir))
                         }, SIMPLIFY = FALSE)
    ## genoimc location to cigar
    MN <- rep(c("M", "N"), readsLen)
    one_zero <- rep(c(1, -1), readsLen)
    cigar <- rep(paste0(readsLen, "M"), length(pos_start))
    dif <- lapply(i_id_range, diff)
    len <- lengths(dif)
    cigar[i_mult_id] <- mapply(dif, len, FUN=function(.ele, l){
      id <- seq.int(l)
      paste(abs(.ele)+one_zero[id], MN[id],
            sep="", collapse = "")
    })
    #njunc <- abs(pos_start$i_id - pos_stop$i_id)
    reads <- GAlignments(seqnames = seqnames(pos_start),
                         strand = strand(pos_start),
                         cigar = cigar,
                         pos = st,
                         names = paste(as.character(seqnames(.sample)),
                                       unlist(sapply(runLength(seqnames(.sample)),
                                              seq.int, simplify = FALSE)),
                                       sep = "_"),
                         tx = as.character(seqnames(.sample)))
    seqinfo(reads) <- seqinfo(txdb)[seqlevels(reads)]

    if(!mis_genome){
      seqinfo(reads) <- seqinfo(genome)[seqlevels(reads)]
      if(includeReadsSeq){
        kgr1 <- kgr <- as(reads, "GRanges")
        strand(kgr) <- "*"
        seqs <- getSeq(genome, kgr)
        seqs <- sequenceLayer(seqs, cigar(reads),
                              from = "reference",
                              to="query-after-soft-clipping")
        seqs[strand(kgr1)=="-"] <- reverseComplement(seqs[strand(kgr1)=="-"])
        mcols(reads)$seq <- seqs
      }
    }else{
      if(any(is.na(seqlengths(reads)))){
        xl <- seqlengths(reads)
        max_reads <- max(end(reads))
        xl[is.na(xl)] <- max_reads
        seqlengths(reads) <- xl
      }
    }
    reads
  })
  ## validation
  # tx <- cds[cds$tx_name %in% "ENSDART00000166393"]
  # re <- start(p[[1]][seqnames(p[[1]]) %in% "ENSDART00000166393"])
  # seq <- getSeq(genome, tx)
  # seq0 <- paste(as.character(seq), collapse="")
  # seq01 <- sapply(re, function(.re) substr(seq0, start = .re, stop = .re+readsLen-1))
  # grn <- gl[[1]]
  # seq2 <- as.character(mcols(grn[mcols(grn)$tx %in% "ENSDART00000166393"])$seq)
  # all(seq01==seq2)
  # tx <- cds[cds$tx_name %in% "ENSDART00000159144"]
  # re <- start(p[[1]][seqnames(p[[1]]) %in% "ENSDART00000159144"]) #ENSDART00000164359
  # seq <- getSeq(genome, tx)
  # seq0 <- paste(as.character(seq), collapse="")
  # seq01 <- sapply(re, function(.re) substr(seq0, start = .re, stop = .re+readsLen-1))
  # seq2 <- as.character(mcols(grn[mcols(grn)$tx %in% "ENSDART00000159144"])$seq)
  # all(seq01==seq2)
  if(!missing(outPath)){
    dir.create(outPath)
    mapply(gl, samples, FUN=function(dat, name){
      tryCatch(
        export(dat, file.path(outPath, paste0(name, ".bam"))),
        error = function(e){
          message("Got error when export the reads.")
          message(e)
        }
      )
    })
  }
  return(invisible(gl))
}
