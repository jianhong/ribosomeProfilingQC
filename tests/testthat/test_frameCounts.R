test_that("frameCounts works not correct", {
  txdb <- makeTxDbFromGFF(system.file("extdata",
                                      "Danio_rerio.GRCz10.91.chr1.gtf.gz",
                                      package="ribosomeProfilingQC"),
                          organism = "Danio rerio",
                          chrominfo = seqinfo(Drerio)["chr1"],
                          taxonomyId = 7955)
  CDS <- prepareCDS(txdb)
  ENSDART00000166393 <- CDS[CDS$tx_name =="ENSDART00000166393"]
  ss <- start(ENSDART00000166393)
  st <- end(ENSDART00000166393)
  a <- unlist(mapply(seq, ss, st, SIMPLIFY = FALSE),
              use.names = FALSE)
  a <- a[(seq_along(a)-1) %% 3 == 0]
  p <- seq(.2, 1, by = .2)
  for(i in p){
    pcs <- GRanges(seqnames = "chr1",
                   ranges = IRanges(sample(a, size=length(a)*i),
                                    width=1),
                   strand = "+")
    pcs <- assignReadingFrame(pcs, CDS)
    fc <- frameCounts(pcs, level = "tx",
                      coverageRate = TRUE, frame0only = TRUE)
    fc <- unname(fc["ENSDART00000166393"])
    expect_equal(fc, i, tolerance=.001)
  }
})
