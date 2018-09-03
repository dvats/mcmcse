\name{multiESS}
\alias{multiESS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effective Sample Size of a multivariate Markov chain as described in Vats et al. (2015).
}
\description{
Calculate the effective sample size of the Markov chain, using the multivariate dependence structure of the process.
}
\usage{
multiESS(x, covmat = NULL, g = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a matrix or data frame of Markov chain output. Number of rows is the Monte Carlo sample size.}
\item{covmat}{optional matrix estimate obtained using \code{mcse.multi} or \code{mcse.initseq}.}
\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}
\item{...}{arguments for \code{mcse.multi} function. Don't use this if a suitable matrix estimate from \code{mcse.multi} or \code{mcse.initseq} is already obtained.}
.}
\details{
Effective sample size is the size of an iid sample with the same variance
  as the current sample. ESS is given by
  \deqn{\mbox{ESS}= n \frac{|\Lambda|^{1/p}}{ |\Sigma|^{1/p},}}{ESS = n |\Lambda|^{1/p}/ |\Sigma|^{1/p},}
  where \eqn{\Lambda}{\Lambda} is the sample covariance matrix for g and \eqn{\Sigma}{\Sigma} is an estimate of the Monte Carlo standard error for g.
}
\value{
  The function returns the estimated effective sample size.
}
\references{
 Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte Carlo, \emph{arXiv preprint arXiv:1512.07713 (2015)}.
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

multiESS(out)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


\seealso{
  \code{\link{minESS}}, which calculates the minimum effective samples required for the problem.

  \code{\link{ess}} which calculates univariate effective sample size using a Markov chain and a function g.
}
