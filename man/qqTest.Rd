\name{qqTest}
\alias{qqTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QQplot for Markov chains 
}
\description{
QQplot for Markov chains using an estimate of the Markov Chain CLT covariance matrix.
}
\usage{
qqTest(mcse.obj)
}
\arguments{
\item{mcse.obj}{the list returned by the \code{mcse.multi} or \code{mcse.initseq} command}
}
\examples{
library(mAr)
p <- 35
n <- 1e4
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
set.seed(100)
foo <- matrix(rnorm(p^2), nrow = p)
foo <- foo \%*\% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)
  
out <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))

mcse.bm <- mcse.multi(x = out)
qqTest(mcse.bm)

mcse.isadj <- mcse.initseq(x = out, adjust = TRUE)
qqTest(mcse.isadj)
}
