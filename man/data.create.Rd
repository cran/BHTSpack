\name{data.create}
\alias{data.create}

\title{
  Create Synthetic Data
}
\description{
 This function generates synthetic compound data.
}
\usage{
data.create(N, nr, nc, M, p, s=NULL, covrow=NULL, covcol=NULL, c=0.0001, mat=FALSE)
}
\arguments{
  \item{N}{
  Number of compounds per plate.
}
  \item{nr}{
  Number of plate rows.
}
  \item{nc}{
  Number of plate columns.
}
  \item{M}{
  Number of plates.
}
  \item{p}{
  Probability of a compound being a hit.
}
  \item{s}{
  Random seed (for reproducibility perposes). Default is \eqn{NULL}.
}
  \item{covrow}{
  Noise plate row-covariance matrix. Default is \eqn{NULL}.
}
  \item{covcol}{
  Noise plate column-covariance matrix. Default is \eqn{NULL}.
}
  \item{c}{
  Constant for scaling plate noise. Default is \eqn{0.0001}.
}
  \item{mat}{
  Specifies a matrix (\eqn{TRUE}) or a vector (\eqn{FALSE}) plate format. Default is \eqn{FALSE}.
}
}

\value{
  This function returns a list consisting of the following elements:
  \item{Z}{A list of matrices (\eqn{mat=TRUE}) or vectors (\eqn{mat=FALSE}) of compounds.}
  \item{B}{A list of compound indicators specifying a hit (\eqn{1}) or a non-hit (\eqn{0}).}
  \item{I}{A list of compound indicators specifying the muxture component (from \eqn{1} to \eqn{K}).}
}

\examples{
 #See package vignette
}
