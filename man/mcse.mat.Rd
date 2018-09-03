\name{mcse.mat}
\alias{mcse.mat}
\title{Apply \code{mcse} to each column of a matrix or data frame of MCMC samples.}
\usage{
  mcse.mat(x, size = "sqroot", g = NULL,
    method = "bm")
}
\arguments{
  \item{x}{a matrix or data frame with each row being a
  draw from the multivariate distribution of interest.}

  \item{size}{the batch size. The default value is
  \dQuote{\code{sqroot}}, which uses the square root of the
  sample size. \dQuote{\code{cuberoot}} will cause the
  function to use the cube root of the sample size. A
  numeric value may be provided if neither
  \dQuote{\code{sqroot}} nor \dQuote{\code{cuberoot}} is
  satisfactory.}

  \item{g}{a function such that \eqn{E(g(x))} is the
  quantity of interest. The default is \code{NULL}, which
  causes the identity function to be used.}

  \item{method}{the method used to compute the standard
  error. This is one of \dQuote{\code{bm}} (batch means,
  the default),\dQuote{\code{bm}} (weighted batch means),
  \dQuote{\code{bm}} (lugsail batch means),
   \dQuote{\code{obm}} (overlapping batch
  means), \dQuote{\code{tukey}} (spectral variance method
  with a Tukey-Hanning window), or \dQuote{\code{bartlett}}
  (spectral variance method with a Bartlett window).}
}
\value{
  \code{mcse.mat} returns a matrix with \code{ncol(x)} rows
  and two columns. The row names of the matrix are the same
  as the column names of \code{x}. The column names of the
  matrix are \dQuote{\code{est}} and \dQuote{\code{se}}.
  The \eqn{j}th row of the matrix contains the result of
  applying \code{mcse} to the \eqn{j}th column of \code{x}.
}
\description{
  Apply \code{mcse} to each column of a matrix or data
  frame of MCMC samples.
}
\seealso{
  \code{\link{mcse}}, which acts on a vector.
  \code{\link{mcse.multi}}, for a multivariate estimate of the Monte Carlo standard error.
  \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which
  compute standard errors for quantiles.
}

