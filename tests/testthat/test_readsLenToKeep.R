test_that("readsLenToKeep or summaryReadsLength works not correct", {
  set.seed(1)
  prob <- c(.01, .01, .05, .1, .77, .05, .01)
  reads <- GRanges("chr1", ranges=IRanges(seq.int(10000), width=1),
                   qwidth=sample(25:31, size = 10000, replace = TRUE,
                                 prob = prob))
  sr <- summaryReadsLength(reads, widthRange = 25:31, plot = FALSE)
  expect_equivalent(as.numeric(sr),
                    sort(prob, decreasing=TRUE),
                    tolerance=.01)
  rtk <- readsLenToKeep(sr, cutoff = .8)
  expect_equivalent(rtk, c(28, 29))
  rtk <- readsLenToKeep(sr, cutoff = .7)
  expect_equivalent(rtk, 29)
})
