% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TPCtools.R
\name{nbcurve}
\alias{nbcurve}
\title{Re-parameterized Norberg functions}
\usage{
nbcurve(x, topt, w, a, b)
}
\arguments{
\item{x}{Temperature}

\item{topt}{Optimum temperature}

\item{w}{Thermal niche width}

\item{a}{Affects the y-intercept}

\item{b}{Exponential scaling}
}
\value{
Predicted exponential growth rate at temperature x
}
\description{
With this formulation, optimum temperature is explicitly the temperature at which
a species achieves its highest growth rate.
}
