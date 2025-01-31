% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funclib_plot.R
\name{plotPASignals}
\alias{plotPASignals}
\title{plotPASignals}
\usage{
plotPASignals(PAdf, resultpath, bsgenome, signals, colorOfBar)
}
\arguments{
\item{PAdf}{A dataframe that contains the PA information of all reads,
Contains at least three columns: PAID, gene, ftr.}

\item{resultpath}{The specified path to save barplot.}

\item{bsgenome}{Reference genomes of the species of interest.}

\item{signals}{The set of PA signal sequences of interest, such as
signals=c("AATAAA","ATTAAA","AAGAAA","AATATA").}

\item{colorOfBar}{The color of barplot, as input in hexadecimal form,
defaults to "#196BE1".}
}
\value{
A barplot.
}
\description{
Draw the frequency distribution diagram of PA signal.
}
\details{
plot the probability distribution of the occurrence of
  user-specified or default PA signal within the first 50bp base range of PA
  site.
}
\examples{
library("BSgenome.Mmusculus.UCSC.mm10")
bsgenome <- BSgenome.Mmusculus.UCSC.mm10
load("./data/output/PAs.RData")
p <- plotPASignals(PAs, "./inst/extdata/figures/", bsgenome = bsgenome)
}
