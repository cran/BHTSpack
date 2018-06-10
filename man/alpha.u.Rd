\name{alpha.u}
\alias{alpha.u}

\title{
 package internal function
}
\description{
 package internal function
}
\usage{
alpha.u(nu, a0, b0, H)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu}{
  Description
}
  \item{a0}{
  Description
}
  \item{b0}{
  Description
}
  \item{H}{
  Description
}
}

\examples{
  M = 5
  H = 10
  a = 10^-6
  b = 10^-6
  nu = lapply(1:M, function(x){rbeta(H, a, b)})
  alpha.u(nu, a, b, H)
}
