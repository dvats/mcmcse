#####################################################################
### Fast calculation of SVE estimators using Herberle et. al (2017)
#####################################################################


mSVEfft <- function (A, b, method = "bartlett")
{
  n <- nrow(A) # A must be centered matrix
  p <- ncol(A)
  w <- as.vector(lag(1:b, n = n, b = b, method = method)) # calculate lags
  w <- c(1, w[1:(n-1)], 0, w[(n-1):1])  # starting steps from FFT paper
  w <- Re(fftw_r2c(w))            
  FF <- matrix(0, ncol = p, nrow = 2*n)  
  FF[1:n,] <- A    
  if(p > 1)  # multivariate
  {
    FF <- mvfftw_r2c (FF)        
    FF <- FF * matrix(w, nrow = 2*n, ncol = p) 
    FF <- mvfftw_c2r(FF) / (2* n ) 
    return ((t(A) %*% FF[1:n, ]) / n )    
  } else if(p == 1)  ##univariate calls
  {
    FF <- fftw_r2c (FF)        
    FF <- FF * matrix(w, nrow = 2*n, ncol = p) 
    FF <- fftw_c2r(FF) / (2* n ) 
    return ((t(A) %*% FF[1:n]) / n )
  }              

}




#####################################################################
### Main function. Estimates the covariance matrix
### Recommend blather = FALSE for users and TRUE for developers
#####################################################################

#' Multivariate Monte Carlo standard errors for expectations
#' 
#' Function returns the estimate of the covariance matrix in the Markov Chain CLT using batch means
#' or spectral variance methods (with different lag windows). The function also returns the Monte
#' Carlo estimate.
#' 
#' 
#' @param x A matrix or data frame of Markov chain output. Number of rows is the Monte
#' Carlo sample size.
#' @param method Any of \dQuote{\code{bm}},\dQuote{\code{obm}},\dQuote{\code{bartlett}},\dQuote{\code{tukey}}. \dQuote{\code{bm}}
#' represents batch means estimator, \dQuote{\code{obm}} represents the overlapping batch means estimator,
#'   and \dQuote{\code{bartlett}} and \dQuote{\code{tukey}} represent the modified-Bartlett window and
#'   the Tukey-Hanning windows for the spectral variance estimators.
#' @param r The lugsail parameters (\code{r}) that converts a lag window into its lugsail
#'   equivalent. Larger values of \code{r} will typically imply less underestimation of \dQuote{\code{cov}},
#'   but higher variability of the estimator. Default is \code{r = 3} and \code{r = 1,2} are
#'   good choices. \code{r > 5} is not recommended.
#' @param size Represents the batch size in \dQuote{\code{bm}} and the truncation point in \dQuote{\code{bartlett}} and
#'  \dQuote{\code{tukey}}. Default is \code{NULL} which implies that an optimal batch size is calculated
#'  using the \code{batchSize} function. Can take character values of \dQuote{\code{sqroot}} and
#'  \dQuote{\code{cuberoot}} or any numeric value between 1 and n/2. \dQuote{\code{sqroot}} means
#'  size is floor(n^(1/2)) and \dQuote{\code{cuberoot}} means size is floor(n^(1/3)).
#' @param g A function that represents features of interest. \code{g} is applied to each row of \code{x} and
#'   thus \code{g} should take a vector input only. If \code{g} is \code{NULL}, \code{g} is set to be identity, which
#'   is estimation of the mean of the target density.
#' @param adjust Defaults to \code{TRUE}. logical for whether the matrix should automatically be adjusted if unstable.
#' @param blather If \code{TRUE}, returns under-the-hood workings of the package.
#' 
#' @return A list is returned with the following components,
#'  \item{cov}{a covariance matrix estimate.}
#'  \item{est}{estimate of g(x).}
#'  \item{nsim}{number of rows of the input x.}
#'  \item{eigen_values}{eigen values of the estimate cov.}
#'  \item{method}{method used to calculate matrix cov.}
#'  \item{size}{value of size used to calculate cov.}
#'  \item{Adjustment_Used}{whether an adjustment was used to calculate cov.}
#' 
#' @usage mcse.multi(x, method = "bm", r = 3, size = NULL, g = NULL,  
#'                   adjust = TRUE, blather = FALSE)
#' @references 
#' Vats, D., Flegal, J. M., and, Jones, G. L Multivariate output analysis for Markov chain Monte Carlo, 
#' \emph{Biometrika}, \bold{106}, 321–-337.
#' 
#' Vats, D., Flegal, J. M., and, Jones, G. L. (2018) Strong Consistency of multivariate spectral variance 
#' estimators for Markov chain Monte Carlo, \emph{Bernoulli}, \bold{24}, 1860–-1909.
#' 
#' @seealso \code{\link{batchSize}}, which computes an optimal batch size. 
#' \code{\link{mcse.initseq}}, which computes an initial sequence estimator.
#' \code{\link{mcse}}, which acts on a vector. 
#' \code{\link{mcse.mat}}, which applies mcse to each column of a matrix or data frame. 
#' \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard
#' errors for quantiles.
#'  
#' @export
#'  
#' @examples 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu <- c(2, 50)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' out <- BVN_Gibbs(n, mu, sigma)
#' 
#' mcse.bm <- mcse.multi(x = out)
#' mcse.tuk <- mcse.multi(x = out, method = "tukey")
#' 
#' # If we are only estimating the mean of the first component,
#' # and the second moment of the second component
#' 
#' g <- function(x) return(c(x[1], x[2]^2))
#' mcse <- mcse.multi(x = out, g = g)

mcse.multi <- function(x, method = "bm", r=3, size = NULL, g = NULL, adjust = TRUE, blather = FALSE)
{ 
  method <- match.arg(method, choices = c("bm", "obm", "bartlett", "tukey", "lug"))
  if(method == "lug")   # not releaved to the public. Inside option for developers
  {
    method <- "bm"
    r <- 3
  }
  # at some point the method used may be different
  # from the method asked. Blather will output this
  method.used <- method

  c <- 0.5  
  if(!is.numeric(r)) stop("r should be numeric")
  if(!is.numeric(c)) stop("c should be numeric")
  
  if(r > 5) warning("We recommend using r <=5. r = 1,2,3 are standard")
  if(r < 1)  {
    warning("r cannot be less than 1. Setting r = 3")
    r <- 3
  }
  if(c > 1) {
    warning("c cannot be greater than 1. Setting c = 0.5")
    c <- 0.5
  }
  
  # making matrix compatible and applying g
  chain <- as.matrix(x)
  if(!is.matrix(chain) && !is.data.frame(chain))
    stop("'x' must be a matrix or data frame.")
  
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
  
  
  ## Setting dimensions on the mcmc output. 
  n <- dim(chain)[1]
  p <- dim(chain)[2]
  
  if(n < (p+1))
    stop("sample size is insufficient for a Markov chain of this dimension")
  
  # Setting batch size based on chosen option
  if(is.null(size))
  {
    b <- batchSize(x = x, method = method)  # optimal
  }
  else if(size == "sqroot")
  {
    b <- floor(sqrt(n))
  } 
  else if(size == "cuberoot") {
    b <- floor(n^(1/3))
  } 
  else {
    if (!is.numeric(size) || size < 1 || size >= n || floor(n/size) <=1) {
      warning("size is either too large, too small, or not a number. Setting 'size' to the optimal
              value")
      size = batchSize(x = x, method = method)
    }
        
    b <- floor(size)
  }

  a <- floor(n/b)
  if(b == 1 && r != 1)
  {
    r <- 1
  }
  ########################## 
  
  
  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)  # this is based on the full output. not n = a*b
  
  ## Setting matrix sizes to avoid dynamic memory 
  sig.mat <- matrix(0, nrow = p, ncol = p)
  
  message <- ""   # will store some info for blather
  
  if(b < (2*r)) 
  {
    r <- 1       
    message <- paste(message, "estimated batch size is low, lugsail not required")
  }
  
  adjust.used <- FALSE
  if(b == 1)  
  {
    init.mat <- var(chain)
    sig.mat <- init.mat
  }
  else  
  {
    ## Batch Means
    if(method == "bm")  
    {
      init.mat <- mbmC(chain, b)
      sig.mat <- init.mat
      if(r>1) {
        sig.mat <- (1/(1-c))*init.mat - (c/(1-c))*mbmC(chain, floor(b/r))
      }
    }
    
    
    ## Overlapping Batch Means
    if(method == "obm") 
    {
      init.mat <- mobmC(chain, b)
      sig.mat <- init.mat
      if(r>1) {
        sig.mat <- (1/(1-c))*init.mat - (c/(1-c))*mobmC(chain, floor(b/r))
      }
    }
    
    
    ## Bartlett or Tukey
    if((method == "bartlett") || (method == "tukey"))
    {
      chain <- scale(chain, center = mu.hat, scale = FALSE)
      init.mat <-  mSVEfft(A = chain, b = b, method = method)
      sig.mat <- init.mat
      if(r>1) 
      {
        sig.mat <- (1/(1-c))*init.mat - (c/(1-c))*mSVEfft(A = chain, b = floor(b/r), method = method)
      }
    }
    adjust.used <- FALSE
    method.used <- paste("Lugsail ", method, " with r = ", r)
    if(r != 1)
    { 
      if(prod(diag(sig.mat) > 0) == 0)  # If diagonals are negative, cannot use larger values of r
      {
        sig.mat <- init.mat
        method.used <- method
        message <- paste(message, paste("Diagonals were negative with r = ", r,". r = 1 was used.", sep = ""), sep = "")
        adjust.used <- TRUE  #whether an adjustment was made
      }
    }
  }
  
  
  
  sig.eigen <- eigen(sig.mat, only.values = TRUE)$values
  
  
  
  if(adjust) # if adjust is FALSE, may output non PD estimator
  {
    if(min(sig.eigen) <= 0)  #needs an adjustment. No need to adjust is not needed
    {
      adjust.used <- TRUE
      warning("Estimated matrix not positive definite. The chain might be highly correlated or very high dimensional. Consider increasing the sample size. Using the default batch means estimator as a substitute.")
      if(method == "bm")  {
        sig.mat = init.mat
      } else  {
        sig.mat = mcse.multi(x, method = "bm", r=1, size = size, g = g, adjust = FALSE, blather = FALSE)$cov
      }
      sig.eigen <- eigen(sig.mat, only.values = TRUE)$values
    }    
  } 
  
  if(blather)
  {
    value = list("cov" = sig.mat, "est" = mu.hat, "nsim" = n, "eigen_values" = sig.eigen,
                 "method" = method.used, "size" = b, "Adjustment_Used" = adjust.used, "message" = message)
  } else {
    value = list("cov" = sig.mat, "est" = mu.hat, "nsim" = n, "eigen_values" = sig.eigen)
  }
  class(value) = "mcmcse"
  value
  
}

#'
#' Check if the class of the object is mcmcse
#' 
#' @param x The object that is checked to belong to the class mcmcse
#' 
#' @return Boolean variable indicating if the input belongs to the class mcmcse
#' 
#' @examples 
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu <- c(2, 50)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#'
#' out <- BVN_Gibbs(n, mu, sigma)
#' is.mcmcse(mcse.multi(out))
#' 
#' @export
#' 
is.mcmcse <- function(x) inherits(x, "mcmcse")






