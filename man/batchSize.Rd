\name{batchSize}
\alias{batchSize}
\docType{package}
\title{
Batch size (truncation point) selection
}
\description{
Function returns the optimal batch size (or truncation point) for a given chain and method. 
}
\usage{
batchSize(x, method = "bm", g = NULL)	
}
\arguments{
\item{x}{a matrix or data frame of Markov chain output. Number of rows is the Monte Carlo sample size.}

\item{method}{any of \code{``bm'', ``obm'', ``bartlett'', ``tukey''}. \code{``bm''} represents batch means estimator, \code{``obm''} represents the overlapping batch means estimator, and \code{``bartlett''} and \code{``tukey''} represent the modified-Bartlett window and the Tukey-Hanning windows for the spectral variance estimators. }

\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}

}
\seealso{
  \code{\link{mcse.multi}}, which calls on \code{batchSize}.
  \code{\link{mcse}}, which calls on \code{batchSize}.
}

\value{
A value of the optimal batch size is returned.
}
\references{
 Liu, Y., Vats, D., and Flegal, J. M. Batch size selection for variance estimators in MCMC, \emph{arXiv preprint arXiv:1804.05975  (2019)}.
}
\examples{
library(mAr)
p <- 3
n <- 1e3
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
set.seed(100)
foo <- matrix(rnorm(p^2), nrow = p)
foo <- foo \%*\% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)
  
out <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))

batchSize(out)
batchSize(out, method = "obm")
batchSize(out, method = "bartlett")
}
