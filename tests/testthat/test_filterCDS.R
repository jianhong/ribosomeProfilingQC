test_that("filterCDS works not correct", {
  txdb <- makeTxDbFromGFF(system.file("extdata",
                                      "Danio_rerio.GRCz10.91.chr1.gtf.gz",
                                      package="ribosomeProfilingQC"),
                          organism = "Danio rerio",
                          chrominfo = seqinfo(Drerio)["chr1"],
                          taxonomyId = 7955)
  CDS <- prepareCDS(txdb)
  cds <- filterCDS(CDS, sizeCutoff = 200)
  cds_width <- cds$wid.cumsum[cds$isLastExonInCDS]
  expect_true(all(cds_width>200))
  cds <- filterCDS(CDS, sizeCutoff = 500)
  cds_width <- cds$wid.cumsum[cds$isLastExonInCDS]
  expect_true(all(cds_width>500))
})
