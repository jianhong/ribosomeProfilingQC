test_that("getPsiteCoordinates works not correct", {
  txdb <- makeTxDbFromGFF(system.file("extdata",
                                      "Danio_rerio.GRCz10.91.chr1.gtf.gz",
                                      package="ribosomeProfilingQC"),
                          organism = "Danio rerio",
                          chrominfo = seqinfo(Drerio)["chr1"],
                          taxonomyId = 7955)

  CDS <- prepareCDS(txdb)
  tmpbam <- tempdir()
  ## test from 5'end
  simulateRPF(txdb, outPath=tmpbam, genome=Drerio,
              samples="s13", readsPerSample = 1e4,
              psite = 13, readsLen = 28,
              frame0 = 1, frame1 = 0, frame2 = 0)
  bamfile <- BamFile(file.path(tmpbam, "s13.bam"))
  pc <- getPsiteCoordinates(bamfile = bamfile,
                            bestpsite = 13,
                            anchor = "5end")
  pc <- assignReadingFrame(pc, CDS=CDS)
  m <- table(pc$readingFrame)
  expect_equal(unname(m["0"]/sum(m)), 1, tolerance=.01)
  ## test from 3'end
  simulateRPF(txdb, outPath=tmpbam, genome=Drerio,
              samples="s14", readsPerSample = 1e4,
              psite = 14, readsLen = 29,
              frame0 = 1, frame1 = 0, frame2 = 0)
  f <- mergeBam(files = file.path(tmpbam, c("s13.bam", "s14.bam")),
                destination = file.path(tmpbam, "m.bam"),
                overwrite = TRUE)
  indexBam(f)
  bamfile <- BamFile(f)
  pc <- getPsiteCoordinates(bamfile = bamfile,
                            bestpsite = 28-13+1,
                            anchor = "3end")
  pc <- assignReadingFrame(pc, CDS=CDS)
  m <- table(pc$readingFrame)
  expect_equal(unname(m["0"]/sum(m)), 1, tolerance=.01)
})
