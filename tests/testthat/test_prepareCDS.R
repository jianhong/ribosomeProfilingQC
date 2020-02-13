test_that("prepareCDS works not correct", {
  txdb <- loadDb(system.file("extdata", "danRer10.chr1.txdb",
                             package="ribosomeProfilingQC"))
  cds <- cds(txdb)
  CDS <- prepareCDS(txdb)
  ol <- findOverlaps(cds, CDS, type = "equal")
  expect_true(all(seq_along(cds) %in% queryHits(ol)))
  expect_true(all(seq_along(CDS) %in% subjectHits(ol)))
  ## check isFirstExonInCDS
  cds <- cdsBy(txdb, by = "tx")
  cds <- unlist(range(cds))
  cds_5 <- promoters(cds, upstream = 0, downstream = 1)
  CDS_isFirst <- promoters(CDS[CDS$isFirstExonInCDS],
                           upstream = 0, downstream = 1)
  ol <- findOverlaps(CDS_isFirst, cds_5, type = "equal")
  expect_true(all(seq_along(CDS_isFirst) %in% queryHits(ol)))
  expect_true(all(seq_along(cds_5) %in% subjectHits(ol)))
  ## check isLastExonInCDS
  switchStrand <- function(.ele){
    levels(strand(.ele)) <- c("-", "+", "*")
    .ele
  }
  end3 <- function(.ele){
    promoters(switchStrand(.ele),
              upstream = 0, downstream = 1)
  }
  cds_3 <- end3(cds)
  CDS_isLast <- end3(CDS[CDS$isLastExonInCDS])
  ol <- findOverlaps(CDS_isLast, cds_3, type = "equal")
  expect_true(all(seq_along(CDS_isLast) %in% queryHits(ol)))
  expect_true(all(seq_along(cds_3) %in% subjectHits(ol)))
})
