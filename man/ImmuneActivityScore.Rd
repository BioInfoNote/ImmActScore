% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImmuneActivityScore.R
\name{ImmuneActivityScore}
\alias{ImmuneActivityScore}
\title{@ ImmuneActivityScore}
\usage{
ImmuneActivityScore(
  expr,
  sampleNumber,
  permTimes,
  type.of.data,
  format.of.file
)
}
\arguments{
\item{expr}{TPM matrix derived from scRNAseq,
or normalized gene expression matrix from bulk samples;
rows are genes, cols are cells or samples.}

\item{sampleNumber}{Numeric value indicating number of samples in profile.}

\item{permTimes}{Numeric value indicating times of permutation, default by 100.}

\item{type.of.data}{Character indicating source of expression data, 'Microarray'or 'RNA-seq'.}

\item{format.of.file}{Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.}
}
\description{
The main function to calculating immune activity levels based on ssGSEA algorithm,
both for multiple samples and single sample profile.
}
\note{
R package 'pheatmap' is required.
}
\author{
Liwen Xu
}
