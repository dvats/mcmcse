###################################################
### estimates the "optimal" batch size using the  
### parametric ethod of Liu et.al.Uses thresholding
### on aR coefficients and option to use only the
### tail of the chain for acf calculation
###################################################

#####################################################################
### Functions estimates the "optimal" batch size using the parametric
### method of Liu et.al
#####################################################################

#' Batch size (truncation point) selection
#' 
#' Function returns the optimal batch size (or truncation point) for a given chain and method.
#' 
#' @details batchSize fits a stationary autoregressive process to approximate the marginals of the Markov chain, selecting the order of
#' the process as the one with the maximum AIC among the models with coefficients greater than a threshold.
#' 
#' @details Final update of batchSize using thresholding on aR coefficients and option to use only the tail of the chain for acf
#' calculation.
#' 
#' 
#' @param x A matrix or data frame of Markov chain output. Number of rows is the Monte
#'   Carlo sample size.
#' @param method Any of \dQuote{\code{bm}},\dQuote{\code{obm}},\dQuote{\code{bartlett}},\dQuote{\code{tukey}}. \dQuote{\code{bm}}
#'   represents batch means estimator, \dQuote{\code{obm}} represents the overlapping batch means estimator,
#'   and \dQuote{\code{bartlett}} and \dQuote{\code{tukey}} represent the modified-Bartlett window and
#'   the Tukey-Hanning windows for the spectral variance estimators.
#' @param g A function that represents features of interest. \code{g} is applied to each row of x and
#'   thus \code{g} should take a vector input only. If \code{g} is \code{NULL}, \code{g} is set to be identity, which
#'   is estimation of the mean of the target density.
#' @param fast Boolean variable for fast estimation using a subset of the Markov chain.
#' acf.
#'   
#' @return A value of the optimal batch size (truncation) is returned.
#' 
#' @references 
#' Liu, Y., Vats, D., and Flegal, J. M. (to appear 2021) Batch size selection for variance estimators in MCMC, Methodology and
#' Computing in Applied Probability .
#' 
#' @seealso \code{\link{mcse.multi}}, which calls on batchSize. \code{\link{mcse}}, which calls on batchSize.
#' 
#' @export
#' 
#' @examples 
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' out = BVN_Gibbs(n, mu, sigma)
#' batchSize(out)
#' batchSize(out, method = "obm")
#' batchSize(out, method = "bartlett")
#' 


batchSize <- function(x, method = c("bm", "obm", "bartlett", "tukey", "sub"), g = NULL, fast = TRUE) {
  method = match.arg(method)

  chain = as.matrix(x)
  if(!is.matrix(chain) && !is.data.frame(chain))
    stop("'x' must be a matrix or data frame.")

  p <- ncol(chain)
  n <- sum(!is.na(chain[,1])) # number of non-missing rows

  if (is.function(g)) 
  {
    chain <- apply(chain, 1, g)
    
    if(is.vector(chain))
    {
      chain <- as.matrix(chain)
    }else
    {
      chain <- t(chain)
    }
  }
  
  if(n < (p+1))
    stop("sample size is insufficient for a Markov chain of this dimension")
  
  if(!is.matrix(chain) && !is.data.frame(chain))
    stop("'x' must be a matrix or data frame.")
  
  order.max <- min(p, n - 1L, floor(10 * log10(n))) # Maximum order up to which AR is fit
  xacf = matrix(, nrow = order.max+1, ncol = p)
  
  if(fast)  {                                       # Use only the tail of the chain to calculate acf
    last = min(n, 5e4)  # lenght of tail used
    chain2 = as.matrix(chain[(n-last+1):n,])
    xacf = sapply(1:p, function(i) acf(chain2[,i], type = "covariance", lag.max = order.max, plot = FALSE,
                                       demean=TRUE, na.action = na.pass)$acf)
  }

  else  {                                           # use the entire chain for acf calculation
    xacf = sapply(1:p, function(i) acf(chain[,i], type = "covariance", lag.max = order.max, plot = FALSE,
                                       demean=TRUE, na.action = na.pass)$acf)
  }

  threshold = qnorm((1.95)/2)/sqrt(n)              # threshold used in confidence interaval calculation
  b = batchsize_cpp(n, p, xacf, order.max, method, threshold)
  b <- min(b, floor(n / (p + 1)))
  if(n > 10)
    b = min(b, floor(n/10))
  return(b)
}
