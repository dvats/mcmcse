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
Version: \tab 1.5-0\cr
Date: \tab 2021-08-29\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
James M. Flegal <jflegal@ucr.edu>,\cr  
John Hughes, \cr   
Dootika Vats, <dootika@iitk.ac.in>,\cr   
Ning Dai, \cr  
Kushagra Gupta, and\cr  
Uttiya Maji  

Maintainer: Dootika Vats <dootika@iitk.ac.in>
}
\references{

  Dai, N and Jones, G.L. (2017)  Multivariate initial sequence estimators in Markov chain Monte Carlo, \emph{ Journal of Multivariate Analysis}, \bold{159}, 184-199.
  
  Flegal, J. M. (2012) Applicability of subsampling
  bootstrap methods in Markov chain Monte Carlo. In
  Wozniakowski, H. and Plaskota, L., editors, \emph{Monte
  Carlo and Quasi-Monte Carlo Methods 2010}, pages 363-372.
  Springer-Verlag.

  Flegal, J. M. and Jones, G. L. (2010) Batch means and
  spectral variance estimators in Markov chain Monte Carlo.
  \emph{The Annals of Statistics}, \bold{38}, 1034-1070.

  Flegal, J. M. and Jones, G. L. (2011) Implementing Markov
  chain Monte Carlo: Estimating with confidence. In Brooks,
  S., Gelman, A., Jones, G. L., and Meng, X., editors,
  \emph{Handbook of Markov Chain Monte Carlo}, pages
  175-197. Chapman & Hall/CRC Press.

 Doss, C. R., Flegal, J. M., Jones, G. L., and Neath, R. C. (2014). Markov chain Monte Carlo estimation of quantiles. \emph{Electronic Journal of Statistics}, \bold{8}, 2448-2478.

  Gong, L., and Flegal, J. M. A practical sequential stopping rule for high-dimensional Markov chain Monte Carlo. \emph{Journal of Computational and Graphical Statistics}, \bold{25},  684-700.
  
  Heberle, J., and Sattarhoff, C. (2017). A fast algorithm for the computation of HAC covariance matrix estimators. \emph{Econometrics}, \bold{5}, 9.
  
  Jones, G. L., Haran, M., Caffo, B. S. and Neath, R.
  (2006) Fixed-width output analysis for Markov chain Monte
  Carlo. \emph{Journal of the American Statistical
  Association}, \bold{101}, 1537-1547.

  Liu, Y., Vats, D., and Flegal, J. M. (2021) Batch size selection for variance estimators in MCMC, \emph{Methodology and Computing in Applied Probability}, to appear.

  Vats, D., Flegal, J. M., and, Jones, G. L Multivariate output analysis for Markov chain Monte Carlo, \emph{Biometrika}, \bold{106}, 321-337.

  Vats, D., Flegal, J. M., and, Jones, G. L. (2018) Strong Consistency of multivariate spectral variance estimators for Markov chain Monte Carlo, \emph{Bernoulli}, \bold{24}, 1860-1909.
  
}

\examples{
n <- 1e3
mu = c(2, 50)
sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
out = BVN_Gibbs(n, mu, sigma)

multiESS(out)
ess(out)
mcse.mat(out)

mcse.bm <- mcse.multi(x = out)
mcse.tuk <- mcse.multi(x = out, method = "tukey")

}
