\name{mcse.initseq}
\alias{mcse.initseq}
\docType{package}
\title{
Multivariate Monte Carlo standard errors for expectations with the initial sequence method of Dai and Jones (2017).
}
\description{
Function returns the estimate of the covariance matrix in the Markov Chain central limit theorem using initial sequence method. This method is designed to give an asymptotically conservative estimate of the Monte Carlo standard error.
}
\usage{
mcse.initseq(x, g = NULL, adjust = FALSE, blather = FALSE)
}

\arguments{
\item{x}{a matrix or data frame of Markov chain output. Number of rows is the Monte Carlo sample size.}

\item{adjust}{logical; if TRUE, an adjustment is made to increase slightly the eigenvalues of the initial sequence estimator. The default is FALSE.}	

\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}

\item{blather}{if TRUE, outputs under  the hood information about the function.}
}

\value{
A list is returned with the following components,
\item{cov}{a covariance matrix estimate using intial sequence method.}

\item{cov.adj}{a covariance matrix estimate using adjusted initial sequence method if the input \code{adjust}=TRUE.}

\item{est}{estimate of g(x).}
\item{nsim}{number of rows of the input \code{x}. Only if blather = TRUE}
\item{adjust}{logical of whether an adjustment was made to the initial sequence estimator. Only if blather = TRUE}
}

\seealso{
\code{initseq{mcmc}}, which is a different univariate initial sequence estimator.
\code{\link{mcse}}, which acts on a vector.
\code{\link{mcse.mat}}, which applies \code{mcse} to each column of a matrix or data frame.
\code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard errors for quantiles.
\code{\link{mcse.multi}}, which estimates the covariance matrix in the Markov Chain CLT using batch means or spectral variance methods.
}
\references{
  Dai, N and Jones, G.L. (2017+)  Multivariate initial sequence estimators in Markov chain Monte Carlo, \emph{ Journal of Multivariate Analysis}.
}
\examples{
library(mAr)
p <- 3
n <- 1000
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
set.seed(100)
foo <- matrix(rnorm(p^2), nrow = p)
foo <- foo \%*\% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)

dat <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))

out.mcse <- mcse.initseq(x = dat)
out.mcse.adj  <- mcse.initseq(x = dat,adjust = TRUE)

# If we are only estimating the mean of the first component, 
# and the second moment of the second component
g <- function(x) return(c(x[1], x[2]^2))
out.g.mcse <- mcse.initseq(x = dat, g = g)
}

