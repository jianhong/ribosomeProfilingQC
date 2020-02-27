txdb <- loadDb(system.file("extdata", "danRer10.chr1.txdb", package="ribosomeProfilingQC"))
CDS <- prepareCDS(txdb)
tmpbam <- tempdir()
bamfilenames <- character(0)
samples <- character(0)
yieldSize <- 10000000

for(psite in c(12, 13, 14)){
  bamfilename <- tempfile(tmpdir = tmpbam, fileext = ".bam")
  sample <- sub(".bam", "", basename(bamfilename))
  reads <- simulateRPF(txdb, tmpbam, samples=sample, genome=Drerio,
                       readsPerSample = 1e5, psite = psite,
                       frame0=.90, frame1=.05, frame2=.05)
  bamfilenames[psite] <- bamfilename
  samples[psite] <- sample
}

bamfile <- lapply(bamfilenames[c(12, 13, 14)], BamFile, yieldSize = yieldSize)
names(bamfile) <- c(12, 13, 14)

test_that("estimatePsite works not correct", {
  for(psite in c(12, 13, 14)){
    estpsite <- estimatePsite(bamfile[[as.character(psite)]], CDS, Drerio)
    bestpsite <- bestPsite(estpsite)
    expect_equal(bestpsite, psite)
  }
})

test_that("readsEndPlot works not correct", {
  for(psite in c(12, 13, 14)){
    h <- readsEndPlot(bamfile[[as.character(psite)]], CDS, toStartCodon=TRUE)
    dist <- as.numeric(names(h))
    gp <- c(dist[dist<0&dist>-10] %% 3, (dist[dist>0]-1) %% 3)
    gph <- rowsum(h[dist>-10], gp)
    m <- c(13, 12, 14)[which.max(gph)]
    expect_equal(m, psite)
  }
})

pcs <- mapply(getPsiteCoordinates, bamfile, c(12, 13, 14))

test_that("summaryReadsLength works not correct", {
  for(psite in c(12, 13, 14)){
    readsLen <- summaryReadsLength(pcs[[as.character(psite)]])
    expect_equal(readsLen, c("28"=1))
  }
})

pcs <- lapply(pcs, assignReadingFrame, CDS=CDS)

test_that("assignReadingFrame works not correct", {
  for(psite in c(12, 13, 14)){
    m <- table(pcs[[as.character(psite)]]$readingFrame)
    expect_gt(m["0"], m["1"])
    expect_gt(m["0"], m["2"])
    m <- plotFrameDensity(pcs[[as.character(psite)]])
    m <- round(m)
    expect_true(abs(m[2]-m[3])<=1)
    expect_true(abs(m[2]-5)<=1)
  }
})

pcs <- lapply(pcs, readsDistribution, txdb=txdb, plot=FALSE)
test_that("readsDistribution works not correct", {
  for(psite in c(12, 13, 14)){
    cs <- mcols(pcs[[as.character(psite)]])[, c("Intron",
                                                "upstream",
                                                "downstream",
                                                "InterGenic")]
    cs <- as.data.frame(cs)
    cs <- colSums(cs)
    expect_true(all(cs==0))
  }
})

