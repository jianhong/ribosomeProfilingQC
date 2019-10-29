test_that("shiftReads works not correct", {
  gal <- GAlignments(seqnames=Rle("chr1", 6), pos=1:6,
    cigar=c("47M", "50M", "50M", "50M", "50M", "50M"),
    strand=Rle(strand(c("+", "-", "+", "-")), c(1, 2, 2, 1)),
    qname=tail(letters, 6),
    isize=c(180, -180, -265, 265, 185, -185))
  shifted <- ribosomeProfilingQC:::shiftReads(gal, shift = 12L)
  shifted2 <- ribosomeProfilingQC:::shiftReads(gal, shift = 11L)
  ns <- as.character(strand(gal))=="+"
  expect_equal(start(shifted)[ns], start(gal)[ns]+12)
  expect_equal(end(shifted)[!ns], end(gal)[!ns]-12)
  expect_equal(width(shifted), width(gal)-12)
  expect_equal(start(shifted2)[ns], start(gal)[ns]+11)
  expect_equal(end(shifted2)[!ns], end(gal)[!ns]-11)
  expect_equal(width(shifted2), width(gal)-11)
})
