% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TPCtools.R
\name{nbcurve.legacy}
\alias{nbcurve.legacy}
\title{Norberg function}
\usage{
nbcurve.legacy(x, copt, w, a, b)
}
\arguments{
\item{x}{Temperature}

\item{copt}{Competitive optimum temperature}

\item{w}{Thermal niche width}

\item{a}{ln(Exponential intercept)}

\item{b}{Exponential scaling}
}
\value{
Predicted exponential growth rate at temperature x
}
\description{
See Norberg et al. 2001, Thomas et al. 2012
}
\details{
Note that with the original formulation, 'Optimum temperature' is the temperature at
which growth rate coincides with the exponential portion of this equation. Only if b
is exactly 0 will this match the temperature at which a species achieves its highest
growth rate. See Thomas et al. 2012 for more on this, and \code{nbcurve()} for an
alternative formulation.
}
