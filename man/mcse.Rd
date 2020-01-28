\name{mcse}
\alias{mcse}
\title{Compute Monte Carlo standard errors for expectations.}
\usage{
  mcse(x, size = NULL, g = NULL, r = 3,
    method = "bm",
    warn = FALSE)
}
\arguments{
  \item{x}{a vector of values from a Markov chain.}

  \item{size}{represents the batch size in ``bm'' and the truncation point in ``bartlett'' and ``tukey''. Default is \code{NULL} which implies that an optimal batch size is calculated using the \code{batchSize()} function. Can take character values of \code{``sqroot''} and \code{``cuberoot''} or any numeric value between 1 and n/2. \code{``sqroot''} means size is floor(n^(1/2)) and ``cuberoot'' means size is floor(n^(1/3)).}

  \item{g}{a function such that \eqn{E(g(x))} is the
  quantity of interest. The default is \code{NULL}, which
  causes the identity function to be used.}

\item{method}{any of \code{``bm'', ``obm'', ``bartlett'', ``tukey''}. \code{``bm''} represents batch means estimator, \code{``obm''} represents overlapping batch means estimator with, \code{``bartlett''} and \code{``tukey''} represents the modified-Bartlett window and the Tukey-Hanning windows for  spectral variance estimators. }

\item{r}{the lugsail parameter that converts a lag window into its lugsail equivalent. Larger values of \code{``r''} will typically imply less underestimation of \code{``cov''},  but higher variability of the estimator. Default is \code{``r = 3''} and \code{``r = 1,2''} are good choices. \code{``r > 5''} is not recommended. Non-integer values are ok.}



  \item{warn}{a logical value indicating whether the
  function should issue a warning if the sample size is too
  small (less than 1,000).}
}
\value{
  \code{mcse} returns a list with two elements:
  \item{est}{an estimate of \eqn{E(g(x))}.} \item{se}{the
  Monte Carlo standard error.}
}
\description{
  Compute Monte Carlo standard errors for expectations.
}
\examples{
# Create 10,000 iterations of an AR(1) Markov chain with rho = 0.9.

n = 10000
x = double(n)
x[1] = 2
for (i in 1:(n - 1))
    x[i + 1] = 0.9 * x[i] + rnorm(1)

# Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using batch means.

mcse(x)
mcse.q(x, 0.1)
mcse.q(x, 0.9)

# Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using overlapping batch means.

mcse(x, method = "obm")
mcse.q(x, 0.1, method = "obm")

# Estimate E(x^2) with MCSE using spectral methods.

g = function(x) { x^2 }
mcse(x, g = g, method = "tukey")
}
\references{
  Flegal, J. M. (2012) Applicability of subsampling
  bootstrap methods in Markov chain Monte Carlo. In
  Wozniakowski, H. and Plaskota, L., editors, \emph{Monte
  Carlo and Quasi-Monte Carlo Methods 2010} (to appear).
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

  Gong, L., and Flegal, J. M. A practical sequential stopping rule for high-dimensional Markov chain Monte Carlo. Journal of Computational and Graphical Statistics (to appear).
  
  Jones, G. L., Haran, M., Caffo, B. S. and Neath, R.
  (2006) Fixed-width output analysis for Markov chain Monte
  Carlo. \emph{Journal of the American Statistical
  Association}, \bold{101}, 1537--1547.

  Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte Carlo, \emph{arXiv preprint arXiv:1512.07713 (2015)}.
}
\seealso{
  \code{\link{mcse.mat}}, which applies \code{mcse} to each
  column of a matrix or data frame.

 \code{\link{mcse.multi}}, for a multivariate estimate of the Monte Carlo standard error.

  \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which
  compute standard errors for quantiles.
}

