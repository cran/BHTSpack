\name{h.pr.u}
\alias{h.pr.u}

\title{
 package internal function
}
\description{
 package internal function
}
\usage{
h.pr.u(z, ih, mu, sigma, pk, K, H, n)
}

\arguments{
  \item{z}{

}
  \item{ih}{

}
  \item{mu}{

}
  \item{sigma}{

}
  \item{pk}{

}
  \item{K}{

}
  \item{H}{

}
  \item{n}{

}
}

\value{
 Describe
}

\examples{
  Nmax = 100
  K = 5
  H = 10
  M = 20
  n = sample(Nmax, M, replace=TRUE)

  z = abs(rnorm(sum(n)))
  mu = abs(rnorm(K))
  sigma = 1/rgamma(n=K, shape=10, rate=10)

  ih = sample(H,sum(n), replace=TRUE)

  pk = lambda.u(rbeta(K, 1, 1))

  h.pr.u(z, ih, mu, sigma, pk, K, H, n)
}
