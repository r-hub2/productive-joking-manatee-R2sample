% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_newtest.R
\name{power_newtest}
\alias{power_newtest}
\title{Power for tests with p values}
\usage{
power_newtest(TS, f, xparam, yparam, TSextra, alpha = 0.05, B = 1000)
}
\arguments{
\item{TS}{routine to calculate test statistics.}

\item{xparam}{values of first parameter under the alternative hypothesis.}

\item{yparam}{values of second parameter under the alternative hypothesis.}

\item{TSextra}{list passed to TS.}

\item{alpha}{=0.05  type I error.}

\item{B}{= 1000 number of simulation runs to estimate the power.}

\item{rxy}{routine that generates data.}
}
\value{
A matrix of power values
}
\description{
This function estimates the power of test routines that calculate p value(s)
}
