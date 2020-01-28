\name{mcse.mat}
\alias{mcse.mat}
\title{Apply \code{mcse} to each column of a matrix or data frame of MCMC samples.}
\usage{
  mcse.mat(x, size = NULL, g = NULL,
    method = "bm")
}
\arguments{
  \item{x}{a matrix or data frame with each row being a
  draw from the multivariate distribution of interest.}

 \item{size}{represents the batch size in ``bm'' and the truncation point in ``bartlett'' and ``tukey''. Default is \code{NULL} which implies that an optimal batch size is calculated using the \code{batchSize()} function. Can take character values of \code{``sqroot''} and \code{``cuberoot''} or any numeric value between 1 and n/2. \code{``sqroot''} means size is floor(n^(1/2)) and ``cuberoot'' means size is floor(n^(1/3)).}

  \item{g}{a function such that \eqn{E(g(x))} is the
  quantity of interest. The default is \code{NULL}, which
  causes the identity function to be used.}


  \item{method}{any of \code{``bm'', ``obm'', ``bartlett'', ``tukey''}. \code{``bm''} represents batch means estimator, \code{``obm''} represents overlapping batch means estimator with, \code{``bartlett''} and \code{``tukey''} represents the modified-Bartlett window and the Tukey-Hanning windows for  spectral variance estimators. }
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

