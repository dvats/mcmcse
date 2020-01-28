\name{mcse.multi}
\alias{mcse.multi}
\docType{package}
\title{
Multivariate Monte Carlo standard errors for expectations.
}
\description{
Function returns the estimate of the covariance matrix in the Markov Chain CLT using batch means or spectral variance methods (with different lag windows). The function also returns the Monte Carlo estimate.
}
\usage{
mcse.multi(x, method = "bm", r = 3, size = NULL, 
	g = NULL, adjust = TRUE, blather = FALSE)	
}
\arguments{
\item{x}{a matrix or data frame of Markov chain output. Number of rows is the Monte Carlo sample size.}

\item{method}{any of \code{``bm'', ``obm'', ``bartlett'', ``tukey''}. \code{``bm''} represents batch means estimator, \code{``obm''} represents overlapping batch means estimator with, \code{``bartlett''} and \code{``tukey''} represents the modified-Bartlett window and the Tukey-Hanning windows for  spectral variance estimators. }

\item{r}{the lugsail parameter that converts a lag window into its lugsail equivalent. Larger values of \code{``r''} will typically imply less underestimation of \code{``cov''},  but higher variability of the estimator. Default is \code{``r = 3''} and \code{``r = 1,2''} are good choices. \code{``r > 5''} is not recommended. Non-integer values are ok.}

\item{size}{represents the batch size in ``bm'' and the truncation point in ``bartlett'' and ``tukey''. Default is \code{NULL} which implies that an optimal batch size is calculated using the \code{batchSize()} function. Can take character values of \code{``sqroot''} and \code{``cuberoot''} or any numeric value between 1 and n/2. \code{``sqroot''} means size is floor(n^(1/2)) and ``cuberoot'' means size is floor(n^(1/3)).}

\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}

\item{adjust}{Defaults to \code{TRUE}. logical for whether the matrix should automatically be adjusted if unstable}

\item{blather}{if \code{TRUE}, returns under-the-hood workings of the package}
}
\seealso{
  \code{\link{batchSize}}, which computes an optimal batch size.
  \code{\link{mcse.initseq}}, which computes an initial sequence estimator.
  \code{\link{mcse}}, which acts on a vector.
  \code{\link{mcse.mat}}, which applies \code{mcse} to each
  column of a matrix or data frame.
  \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which
  compute standard errors for quantiles.
}

\value{
A list is returned with the following components,
\item{cov}{a covariance matrix estimate.}
\item{est}{estimate of g(x).}
\item{nsim}{number of rows of the input \code{x}.}
\item{method}{method used to calculate matrix \code{cov}.}
\item{size}{value of size used to calculate \code{cov}.}
\item{adjust.used}{whether an adjustment was used to calculate \code{cov}.}
}
\references{
 Vats, D., Flegal, J. M., and, Jones, G. L (2019) Multivariate Output Analysis for Markov chain Monte Carlo, \emph{Biometrika}.

  Vats, D., Flegal, J. M., and, Jones, G. L. (2018) Strong Consistency of multivariate spectral variance estimators for Markov chain Monte Carlo, \emph{Bernoulli}.
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

mcse.bm <- mcse.multi(x = out)
mcse.tuk <- mcse.multi(x = out, method = "tukey")

# If we are only estimating the mean of the first component, 
# and the second moment of the second component

g <- function(x) return(c(x[1], x[2]^2))
mcse <- mcse.multi(x = out, g = g)
}
