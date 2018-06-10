\name{ptrace}
\alias{ptrace}

\title{
 Trace (ACF) Plots
}
\description{
 This function outputs trace plots of certain latent variables.
}
\usage{
ptrace(res, var, ndisc, nr, nc, type="trace")
}
\arguments{
  \item{res}{
  An output object from \eqn{bhts()}.
}
  \item{var}{
  Variable for which to display convergence diagnostic plots. Current options are \eqn{mu0} (displaying \eqn{\mu_{01},\ldots,\mu_{0K}}), \eqn{mu1} (displaying \eqn{\mu_{11},\ldots,\mu_{1K}}), \eqn{sigma0} (displaying \eqn{\sigma_{01}^2,\ldots,\sigma_{0K}^2}), \eqn{sigma1} (displaying \eqn{\sigma_{11}^2,\ldots,\sigma_{1K}^2}), \eqn{pk0} (displaying \eqn{\lambda_{1}^{(0)},\ldots,\lambda_{K}^{(0)}}) and \eqn{pk1} (displaying \eqn{\lambda_{1}^{(1)},\ldots,\lambda_{K}^{(1)}}).
}
  \item{ndisc}{
  Number of iterations for which to discard samples.
}
  \item{nr}{
  Number of rows in the resulting composite plot.
}
  \item{nc}{
  Number of columns in the resulting composite plot.
}
  \item{type}{
  Type of convergence diagnostic. Currently implemented are trace plots (default \eqn{type="trace"}) and ACF plots (\eqn{type="acf"})
}
}

\examples{
 #See package vignette
}
