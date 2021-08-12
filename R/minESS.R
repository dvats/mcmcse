#' Minimum effective sample size required for stable estimation as described in Vats et al. (2015)
#' 
#'
#' The function calculates the minimum effective sample size required for a specified relative
#' tolerance level. This function can also calculate the relative precision in estimation for a given
#' estimated effective sample size.
#' 
#' @usage minESS(p, alpha = .05, eps = .05, ess = NULL)
#' 
#' @param x Dimension of the estimation problem.
#' @param alpha Confidence level.
#' @param eps Tolerance level. The eps value is ignored is `ess` is not `NULL`.
#' @param ess Estimated effective sample size. Usually the output value from `multiESS`.
#' 
#' @details 
#' The minimum effective samples required when estimating a vector of length `p`, with \eqn{100(
#' 1-\alpha)\%} confidence and tolerance of \eqn{\epsilon} is \deqn{mESS \geq \frac{2^{2/p} \pi}{(p
#' \Gamma(p/2))^{2/p}} \frac{\chi^{2}_{1-\alpha,p}}{\epsilon^{2}}}.
#' The above equality can also be used to get \eqn{\epsilon} from an already obtained estimate of
#' mESS.
#' 
#' @return By default function returns the minimum effective sample required for a given eps
#' tolerance. If `ess`is specified, then the value returned is the `eps` corresponding to that `ess`.
#' 
#' @references 
#' Gong, L., and Flegal, J. M. A practical sequential stopping rule for high-dimensional Markov chain
#' Monte Carlo. Journal of Computational and Graphical Statistics (to appear).
#' Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte
#' Carlo, arXiv preprint arXiv:1512.07713 (2015).
#' 
#' @seealso \code{\link{multiESS}}, which calculates multivariate effective sample size using a
#' Markov chain and a function g.
#' \code{\link{ess}} which calculates univariate effective sample size using a Markov chain and a
#' function g.
#' 
#' @examples 
#' minESS(p = 5)
#' 
#' @export
#' 

minESS <- function(p, alpha = .05, eps = .05, ess = NULL)
{
  crit <- qchisq(1-alpha, p)
  foo <- 2/p
  
  if(is.null(ess) == TRUE)
  {
    logminESS <- foo*log(2) + log(pi) - foo*log(p) - foo*lgamma(p/2) - 2*log(eps) + log(crit)
    names(logminESS) <- "minESS"
    return(round(exp(logminESS)))
  }else{
    if(is.numeric(ess) == FALSE)
    {
      stop("Only numeric entry allowed for ess")
    }
    logEPS <- .5*foo*log(2) + .5*log(pi) - .5*foo*log(p) - .5*foo*lgamma(p/2) - .5*log(ess) + .5*log(crit)
    names(logEPS) <- "Epsilon"
    return(exp(logEPS))
  }
}