% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readsDistribution.R
\name{readsDistribution}
\alias{readsDistribution}
\title{Plot reads distribution in genomic elements}
\usage{
readsDistribution(
  reads,
  txdb,
  upstreamRegion = 3000,
  downstreamRegion = 3000,
  plot = TRUE,
  ...
)
}
\arguments{
\item{reads}{Output of \link{getPsiteCoordinates}}

\item{txdb}{A TxDb object}

\item{upstreamRegion, downstreamRegion}{The range for promoter region and
downstream region.}

\item{plot}{Plot the distribution or not}

\item{...}{Not use.}
}
\value{
The reads with distribution assignment
}
\description{
Plot the percentage of reads in CDS, 5'UTR, 3'UTR, introns, and
 other elements.
}
\examples{
library(Rsamtools)
bamfilename <- system.file("extdata", "RPF.WT.1.bam",
                           package="ribosomeProfilingQC")
yieldSize <- 10000000
bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
pc <- getPsiteCoordinates(bamfile, bestpsite=11)
pc.sub <- pc[pc$qwidth \%in\% c(29, 30)]
library(GenomicFeatures)
library(BSgenome.Drerio.UCSC.danRer10)
txdb <- makeTxDbFromGFF(system.file("extdata",
          "Danio_rerio.GRCz10.91.chr1.gtf.gz",
          package="ribosomeProfilingQC"),
          organism = "Danio rerio",
          chrominfo = seqinfo(Drerio)["chr1"],
          taxonomyId = 7955)
pc.sub <- readsDistribution(pc.sub, txdb, las=2)
}
