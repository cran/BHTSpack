\name{bhts}
\alias{bhts}

\title{
 Bayesian High-Throughput Screening
}
\description{
 This is the package main function.
}
\usage{
bhts(Z, iters, H, K, mu00=NULL, mu10=NULL, a.alpha, b.alpha, a.tau, b.tau, 
pnorm=FALSE, s=NULL, store=FALSE)
}

\arguments{
  \item{Z}{
  A list of compounds.
}
  \item{iters}{
  Number of iterations to perform.
}
  \item{H}{
  Number of local DP components.
}
  \item{K}{
  Number of global DP components.
}
  \item{mu00}{
  Activity level (mean) of non-hit compounds
}
  \item{mu10}{
  Activity level (mean) of hit compounds
}
  \item{a.alpha}{
  Gamma shape parameter specifying local DP concentration prior.
}
  \item{b.alpha}{
  Gamma rate parameter specifying local DP concentration prior.
}
  \item{a.tau}{
  Gamma shape parameter specifying global DP concentration prior.
}
  \item{b.tau}{
  Gamma rate parameter specifying global DP concentration prior.
}
  \item{pnorm}{
  Plate normalization. If \eqn{TRUE}, each plate is normalized 
  to zero mean and unit variance, prior to analysis. Default is \eqn{FALSE}.
}
  \item{s}{
  Random seed (for reproducibility purposes). Default is \eqn{NULL}.
}
  \item{store}{
  If \eqn{TRUE}, all samples of certain latent variables are stored in the output object. 
  Default is \eqn{FALSE}.
}
}

\value{
  This function returns a list consisting of the following elements:
  \item{hatpai}{A list of vectors of posterior probabilities, 
estimating the probability of a compound being a hit.}
  \item{dat.store}{If \eqn{store=TRUE} (default is \eqn{FALSE}), 
the output contains a list of \eqn{iters}\eqn{\times}\eqn{K} matrices 
of samples. Each matrix contains the samples of a separate latent variable. 
At each iteration, the following six variables are stored in a different row
 of their corrpsponding matrix, \eqn{(\lambda_{1}^{(0)},\ldots,\lambda_{K}^{(0)})}, 
\eqn{(\lambda_{1}^{(1)},\ldots,\lambda_{K}^{(1)})}, \eqn{(\mu_{01},\ldots,\mu_{0K})}, 
\eqn{(\mu_{11},\ldots,\mu_{1K})}, \eqn{(\sigma_{01}^2,\ldots,\sigma_{0K}^2)} 
and \eqn{(\sigma_{11}^2,\ldots,\sigma_{1K}^2)}.}
}

\examples{
  set.seed(1234)
  Nmax = 100
  M = 100
  n = sample(Nmax, M, replace=TRUE)
  Z = lapply(n, function(x){abs(rnorm(x))})
  bhts(Z, iters=100, H=10, K=5, mu00=0, mu10=10, a.alpha=10, b.alpha=5, a.tau=10, b.tau=5)
}
