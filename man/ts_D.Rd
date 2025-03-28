% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ts_D}
\alias{ts_D}
\title{This function calculates the test statistics for discrete data}
\usage{
ts_D(typeTS, x, y, vals, TS, TSextra, adw)
}
\arguments{
\item{typeTS}{format of TS}

\item{x}{discrete data set (counts)}

\item{y}{discrete data set (counts)}

\item{vals}{values of discrete RV}

\item{TS}{routine}

\item{TSextra}{list passed to TS function}

\item{adw}{vector of weights for Anderson-Darling test}
}
\value{
A vector of numbers
}
\description{
This function calculates the test statistics for discrete data
}
\keyword{internal}
