\name{mcmcse-package}
\alias{mcmcse-package}
\alias{mcmcse}
\docType{package}
\title{
Monte Carlo Standard Errors for MCMC
}
\description{
Provides tools for computing Monte Carlo standard
errors (MCSE) in Markov chain Monte Carlo (MCMC) settings. MCSE
computation for expectation and quantile estimators is
supported. The package also provides functions for computing
effective sample size and for plotting Monte Carlo estimates
versus sample size.
}
\details{
\tabular{ll}{
Package: \tab mcmcse\cr
Type: \tab Package\cr
Version: \tab 1.3-2\cr
Date: \tab 2017-07-03\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
James M. Flegal <jflegal@ucr.edu>,\\
John Hughes <j.hughes@ucdenver.edu> \\
Dootika Vats <dootika@iitk.ac.in> and\\
Ning Dai <daixx224@umn.edu>

Maintainer: Dootika Vats <dootika@iitk.ac.in>
}
\references{

  Dai, N and Jones, G.L. (2017+)  Multivariate initial sequence estimators in Markov chain Monte Carlo, \emph{ Journal of Multivariate Analysis}.
  
  Flegal, J. M. (2012) Applicability of subsampling
  bootstrap methods in Markov chain Monte Carlo. In
  Wozniakowski, H. and Plaskota, L., editors, \emph{Monte
  Carlo and Quasi-Monte Carlo Methods 2010}, pages 363--372.
  Springer-Verlag.

  Flegal, J. M. and Jones, G. L. (2010) Batch means and
  spectral variance estimators in Markov chain Monte Carlo.
  \emph{The Annals of Statistics}, \bold{38}, 1034--1070.

  Flegal, J. M. and Jones, G. L. (2011) Implementing Markov
  chain Monte Carlo: Estimating with confidence. In Brooks,
  S., Gelman, A., Jones, G. L., and Meng, X., editors,
  \emph{Handbook of Markov Chain Monte Carlo}, pages
  175--197. Chapman & Hall/CRC Press.

  Flegal, J. M., Jones, G. L., and Neath, R. (2012) Markov
  chain Monte Carlo estimation of quantiles.
  \emph{University of California, Riverside, Technical
  Report}.

  Gong, L., and Flegal, J. M. A practical sequential stopping rule for high-dimensional Markov chain Monte Carlo. \emph{Journal of Computational and Graphical Statistics}, \bold{25},  684–-700.
  
  Jones, G. L., Haran, M., Caffo, B. S. and Neath, R.
  (2006) Fixed-width output analysis for Markov chain Monte
  Carlo. \emph{Journal of the American Statistical
  Association}, \bold{101}, 1537--1547.

  Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte Carlo, \emph{Biometrika}, \bold{106}, 321–-337.

  Vats, D., Flegal, J. M., and, Jones, G. L. (2018) Strong Consistency of multivariate spectral variance estimators for Markov chain Monte Carlo, \emph{Bernoulli}, \bold{24}, 1860–-1909.

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

mcse(out[,1], method = "bart")

mcse.bm <- mcse.multi(x = out)
mcse.tuk <- mcse.multi(x = out, method = "tukey")

}
