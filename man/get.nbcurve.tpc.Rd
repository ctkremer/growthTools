% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TPCtools.R
\name{get.nbcurve.tpc}
\alias{get.nbcurve.tpc}
\title{Fit Norberg curve to growth rate vs. temperature data}
\usage{
get.nbcurve.tpc(
  temperature,
  mu,
  method = "grid.mle2",
  plotQ = FALSE,
  conf.bandQ = TRUE,
  fpath = NA,
  id = NA,
  suppress.grid.mle2.warnings = TRUE,
  ...
)
}
\arguments{
\item{temperature}{Temperature}

\item{mu}{Exponential growth rate}

\item{method}{Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'}

\item{plotQ}{Should regression be visualized?}

\item{conf.bandQ}{Should we calculate a confidence band around the regression? logical.}

\item{fpath}{If visual requested, and valid file path provided here, plot will be saved as a .pdf file. Default is NA.}

\item{id}{Character string providing any information ID'ing the specifc curve being fit; used to label plots, if any are requested. Default is NA.}

\item{suppress.grid.mle2.warnings}{logical; should warnings arising from grid.mle2 invocation be suppressed (TRUE), or displayed (FALSE)? Default is TRUE.}

\item{...}{Additional arguments passed to grid.mle2 (e.g., control=list(maxit=2000))}
}
\description{
Note: this function fits a re-parameterized version of the original Norberg curve (Norberg et al. 2001, Thomas et al. 2012), altered to depend directly on a parameter that provides the true optimum temperature (ie, the temperature at which growth rate is maximal). See \code{nbcurve()} for more.
}
\examples{
data("example_TPC_data")
sp1 <- example_TPC_data \%>\% filter(isolate.id=='CH30_4_RI_03' & dilution==1)
nbcurve.traits<-get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',
plotQ=TRUE,conf.bandQ = TRUE,fpath=NA)
nbcurve.traits

}
