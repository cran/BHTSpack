\name{mu.k.u}
\alias{mu.k.u}

\title{
 package internal function
}
\description{
 package internal function
}
\usage{
mu.k.u(k, ik, z, sigma, mu0)
}

\arguments{
  \item{k}{

}
  \item{ik}{

}
  \item{z}{

}
  \item{sigma}{

}
  \item{mu0}{

}
}

\value{
 Describe
}

\examples{
  K = 5
  n = 100
  z = abs(rnorm(n))
  sigma = sapply(1:K, function(x){1/rgamma(n=1, shape=10, rate=10)})
  mu0 = 0
  ik = sample(K, n, replace=TRUE)
  sapply(1:K, mu.k.u, ik, z, sigma, mu0)
}
