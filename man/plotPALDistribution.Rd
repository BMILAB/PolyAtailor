% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funclib_plot.R
\name{plotPALDistribution}
\alias{plotPALDistribution}
\title{plotPALDistribution}
\usage{
plotPALDistribution(taildf, resultpath, dType, medianPAL)
}
\arguments{
\item{taildf}{A dataframe that contains the tail information of all reads,
Contains at least three columns: PAL, gene, umi, and sample.}

\item{resultpath}{The specified path to save barplot.}

\item{dType}{The tail length affiliation relationship that you want to draw,
which includes three options: "global", "gene", "umi". "global" A stands
for uniform treatment of all tail lengths. "gene" stands for tail length
plotted according to genetic statistics. "umi" stands for tail length
plotted according to umi statistics}

\item{medianPAL}{A Boolean value. T means Annotating peak. F means the
opposite.}
}
\value{
A density curve with peak annotations.
}
\description{
Plot the PAL distribution of gene, umi or global.
}
\details{
Using the PoltPAL family of functions, the tail length distribution
  of the gene or UMI of different sample will be plotted, as well as the
  global tail length distribution.
}
\examples{
load("./data/output/AnnotatedTails.RData")
p1 <- plotPALDistribution(AnnotedTails, "./inst/extdata/figures/", "global", medianPAL = T)
p2 <- plotPALDistribution(AnnotedTails, "./inst/extdata/figures/", "gene", medianPAL = T)
}
