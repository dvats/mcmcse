\name{confRegion}
\alias{confRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence regions (ellipses) for Monte Carlo estimates
}
\description{
Constructs confidence regions (ellipses) from the Markov chain output for the features of interest. Function uses the ellipse package.
}
\usage{
confRegion(mcse.obj, which = c(1,2), level = .95)
}

\arguments{
\item{mcse.obj}{the list returned by the \code{mcse.multi} or \code{mcse.initseq} command}
\item{which}{integer vector of length 2 indicating the component for which to make the confidence ellipse. Chooses the first two by default.}

\item{level}{confidence level for the ellipse}
}
\details{
Returns a matrix of x and y coordinates for the ellipse. Use plot function on the matrix to plot the ellipse
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
mcerror <- mcse.multi(out)

## Plotting the ellipse
plot(confRegion(mcerror), type = 'l')
}

