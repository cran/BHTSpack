\name{sig.k.u}
\alias{sig.k.u}

\title{
 package internal function
}
\description{
 package internal function
}
\usage{
sig.k.u(k, ik, z, mu0, a0, b0)
}

\arguments{
  \item{k}{

}
  \item{ik}{

}
  \item{z}{

}
  \item{mu0}{

}
  \item{a0}{

}
  \item{b0}{

}
}

\value{
 Describe
}

\examples{
  K = 5
  n = 100
  z = abs(rnorm(n))
  mu0 = 0
  ik = sample(K, n, replace=TRUE)
  a0 = 5
  b0 = 5
  sapply(1:K, sig.k.u, ik, z, mu0, a0, b0)
}
