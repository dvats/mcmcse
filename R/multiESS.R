#' Effective Sample Size of a multivariate Markov chain as described in Vats et al. (2015).
#' 
#' Calculate the effective sample size of the Markov chain, using the multivariate dependence
#' structure of the process.
#' 
#' @usage multiESS(x, covmat = NULL, g = NULL, ...).
#' 
#' @param x a matrix or data frame of Markov chain output. Number of rows is the Monte Carlo sample
#' size.
#' @param covmat optional matrix estimate obtained using mcse.multi or `mcse.initseq`.
#' @param g a function that represents features of interest. g is applied to each row of x and
#' thus g should take a vector input only. If g is `NULL`, g is set to be identity, which
#' is estimation of the mean of the target density.
#' @param ... arguments for `mcse.multi` function. Don't use this if a suitable matrix estimate
#' from `mcse.multi` or `mcse.initseq` is already obtained.
#' 
#' @details 
#' Effective sample size is the size of an iid sample with the same variance as the current sample.
#' ESS is given by \deqn{ESS = n\frac{|\Lambda|^{1/p}}{|\Sigma|^{1/p}},} where \eqn{\Lambda} is the
#' sample covariance matrix for g and \eqn{\Sigma} is an estimate of the Monte Carlo standard
#' error for g.
#' 
#' @return The function returns the estimated effective sample size.
#' 
#' @references 
#' Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte
#' Carlo, arXiv preprint arXiv:1512.07713 (2015).
#' 
#' @seealso \code{\link{minESS}}, which calculates the minimum effective samples required for the
#' problem.
#' \code{\link{ess}} which calculates univariate effective sample size using a Markov chain and a
#' function g.
#' 
#' @examples 
#' library(mAr)
#' p <- 3
#' n <- 1e3
#' omega <- 5*diag(1,p)
#' ## Making correlation matrix var(1) model
#' set.seed(100)
#' foo <- matrix(rnorm(p^2), nrow = p)
#' foo <- foo %*% t(foo)
#' phi <- foo / (max(eigen(foo)$values) + 1)
#' out <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))
#' multiESS(out)
#' 
#' library(mvtnorm)
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' p <- 2
#' n <- 1e3
#' mu1 <- 2 
#' mu2 <- 50
#' A <- 1
#' B <- 1
#' rho <- 0.5
#' sigma = matrix(c(A, rho, rho, B), nrow = 2)
#' init = rmvnorm(1, mean = c(mu1, mu2), sigma = sigma) ## Starting from stationarity
#' X <- matrix(0, nrow = n, ncol = p)
#' X[1, ] = init
#' ## Gibbs sampler to generate the Markov chain
#' for (i in 2:n) {
#'  X[i, 1] = rnorm(1, mu1 + (rho / b) * (X[i - 1, 2] - mu2), sqrt(a - (rho ^ 2) / b))
#'  X[i, 2] = rnorm(1, mu2 + (rho / a) * (X[i, 1] - mu1), sqrt(b - (rho ^ 2) / a))
#' }
#' multiESS(X)
#' 
#' @export
#' 

multiESS <- function(x, covmat = NULL, g = NULL, ...)
{

	chain <- as.matrix(x)
	if(!is.matrix(x) && !is.data.frame(x))
	  stop("'x' must be a matrix or data frame.")

	if (is.function(g)) 
	{
	  chain <- apply(x, 1, g)

	  if(is.vector(chain))
	  {
	    chain <- as.matrix(chain)
	  }else
	  {
	    chain <- t(chain)
	  }
	}
	## Setting dimensions on the mcmc output. 
	n = dim(chain)[1]
	p = dim(chain)[2]

	if(is.matrix(covmat))
	{
		var_mat <- cov(chain)
		det.var.p <- prod(eigen(var_mat, only.values = TRUE)$values^(1/p))
		det.covmat.p <- prod(eigen(covmat, only.values = TRUE)$values^(1/p))
		ess <- n*(det.var.p/det.covmat.p)
	} else
	{
		covmat <- mcse.multi(chain, ...)$cov
		var_mat <- cov(chain)

		det.var.p <- prod(eigen(var_mat, only.values = TRUE)$values^(1/p))
		det.covmat.p <- prod(eigen(covmat, only.values = TRUE)$values^(1/p))
		ess <- n*(det.var.p/det.covmat.p)

	}
return(ess)

}