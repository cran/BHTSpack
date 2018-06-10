\name{b.u}
\alias{b.u}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 package internal function
}
\description{
 package internal function
}
\usage{
b.u(hatpai)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hatpai}{
  Description
}
}

\examples{
pai = 0.5
M = 10
H = 10
K = 5
n = 100

z = abs(rnorm(n))

sigma1 = abs(rnorm(K))
sigma0 = abs(rnorm(K))

mu1 = abs(rnorm(K))
mu0 = abs(rnorm(K))

hk0 = matrix(sample(K, M*H, replace=TRUE), M, H)
hk1 = matrix(sample(K, M*H, replace=TRUE), M, H)

nu.h0 = lapply(1:H, function(x){rbeta(1,5,5)})
nu.h1 = lapply(1:H, function(x){rbeta(1,5,5)})

ph0 = lapply(nu.h0, lambda.u)
ph1 = lapply(nu.h1, lambda.u)

b.u(hatpai.u(z, hk1, hk0, ph1, ph0, sigma1, sigma0, mu1, mu0, pai, H, n))
}
