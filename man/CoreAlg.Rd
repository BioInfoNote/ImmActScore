% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImmuneInfiltration.R
\name{CoreAlg}
\alias{CoreAlg}
\title{Core algorithm}
\usage{
CoreAlg(X, y)
}
\arguments{
\item{X}{cell-specific gene expression.}

\item{y}{mixed expression per sample.}
}
\description{
It's a subfunction of CIBERSORT.
SVR algorithm, the core algorithm of CIBERSORT, referenced by CIBERSORT to compute infiltration of immune cells and it can't run alone.
}
\details{
Download this file from https://cibersort.stanford.edu/
CIBERSORT R script v1.03 (last updated 07-10-2015)
Primary Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
There is no change compared to the source code from CIBERSORT.
License: http://cibersort.stanford.edu/CIBERSORT_License.txt
}
