#' Effective Sample Size of a multivariate Markov chain as described in Vats et al. (2015)
#' 
#' Calculate the effective sample size of the Markov chain, using the multivariate dependence
#' structure of the process.
#' 
#' 
#' @param x a matrix or data frame of Markov chain output. Number of rows is the Monte Carlo sample
#' size.
#' @param covmat optional matrix estimate obtained using \code{mcse.multi} or \code{mcse.initseq}.
#' @param g a function that represents features of interest. \code{g} is applied to each row of \code{x} and
#' thus \code{g} should take a vector input only. If \code{g} is \code{NULL}, \code{g} is set to be identity, which
#' is estimation of the mean of the target density.
#' @param ... arguments for \code{mcse.multi} function. Don't use this if a suitable matrix estimate
#' from \code{mcse.multi} or \code{mcse.initseq} is already obtained.
#' 
#' @details 
#' Effective sample size is the size of an iid sample with the same variance as the current sample.
#' ESS is given by \deqn{ESS = n\frac{|\Lambda|^{1/p}}{|\Sigma|^{1/p}},} where \eqn{\Lambda} is the
#' sample covariance matrix for \code{g} and \eqn{\Sigma} is an estimate of the Monte Carlo standard
#' error for \code{g}.
#' 
#' @return The function returns the estimated effective sample size.
#' 
#' @references 
#' Vats, D., Flegal, J. M., and, Jones, G. L Multivariate output analysis for Markov chain Monte Carlo, 
#' \emph{Biometrika}, \bold{106}, 321-337.
#' 
#' @seealso \code{\link{minESS}}, which calculates the minimum effective samples required for the
#' problem.
#' \code{\link{ess}} which calculates univariate effective sample size using a Markov chain and a
#' function g.
#' 
#' @examples 
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu <- c(2, 50)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' out <- BVN_Gibbs(n, mu, sigma)
#'
#' multiESS(out)
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
    if(is.mcmcse(covmat)) {
      eigs_cov = covmat$eigen_values
    } else  {
      eigs_cov = eigen(covmat, only.values = TRUE)$values
    }
  } else
  {
    covmat <- mcse.multi(x, ...)$cov
    eigs_cov = eigen(covmat, only.values = TRUE)$values
  }
  var_mat <- cov(chain)
  
  log.det.var.p <- sum(log(eigen(var_mat, symmetric = TRUE, only.values = TRUE)$values))
  log.det.covmat.p <- sum(log(eigs_cov))
  ess <- n*exp((log.det.var.p - log.det.covmat.p)/p)
  return(ess)
  
}
