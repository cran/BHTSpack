\name{r.fdr}
\alias{r.fdr}

\title{
 Significant Hits
}
\description{
 This function determines significant hits, based on a specified expected FDR.
}
\usage{
r.fdr(res, fdr=0.05)
}
\arguments{
  \item{res}{
  An output object from \eqn{bhts()}.
}
  \item{fdr}{
  Expected FDR (default is 0.05).
}
}

\value{
  This function returns a list consisting of the following elements:
  \item{res}{A data frame containing significant hits and their probabilities.}
  \item{r}{The computed significant hit probability threshold.}
}

\examples{
 #See package vignette
}
