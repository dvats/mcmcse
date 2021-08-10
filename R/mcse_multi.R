# library(Rcpp)
# library(RcppArmadillo)
# sourceCpp("mobmc.cpp")
# sourceCpp("mbmc.cpp")
# sourceCpp("msvec.cpp")

#####################################################################
### Function adjusts a non-positive definite estimator to be PD
#####################################################################

adjust_matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 9/10)
{
  mat.adj <- mat
  adj <- epsilon*N^(-b)
  vars <- diag(mat)
  corr <- cov2cor(mat)
  eig <- eigen(corr)
  adj.eigs <- pmax(eig$values, adj)
  mat.adj <- diag(vars^(1/2))%*% eig$vectors %*% diag(adj.eigs) %*% t(eig$vectors) %*% diag(vars^(1/2))
  return(mat.adj)
}

#####################################################################
### Function approximates each column of the output by an AR(p)
### and estimates Gamma and Sigma to calculate batch size
#####################################################################

arp_approx <- function(x)
{
  
  # Fitting a univariate AR(m) model
  ar.fit <- ar(x, aic = TRUE)

  # estimated autocovariances
  gammas <- as.numeric(acf(x, type = "covariance", lag.max = ar.fit$order, plot = FALSE)$acf)
  spec <- ar.fit$var.pred/(1-sum(ar.fit$ar))^2  #asym variance
  
  if(ar.fit$order != 0)
  {
  foo <- 0
  for(i in 1:ar.fit$order)
  {
    for(k in 1:i)
    {
      foo <- foo + ar.fit$ar[i]*k*gammas[abs(k-i)+1]
    }
  }
  Gamma <- 2*(foo + (spec - gammas[1])/2 *sum(1:ar.fit$order * ar.fit$ar)  )/(1-sum(ar.fit$ar))
  } else{
    Gamma <- 0
  }
  rtn <- cbind(Gamma, spec)
  colnames(rtn) <- c("Gamma", "Sigma")
  return(rtn)
}

#####################################################################
#####################################################################
### Functions below will be exported
#####################################################################
#####################################################################


#####################################################################
### Functions estimates the "optimal" batch size using the parametric
### method of Liu et.al
#####################################################################

#' Batch size (truncation point) selection
#' 
#' Function returns the optimal batch size (or truncation point) for a given chain and method.
#' 
#' @usage batchSize(x, method = "bm", g = NULL)
#' 
#' @param x A matrix or data frame of Markov chain output. Number of rows is the Monte
#'   Carlo sample size.
#' @param method Any of `bm`,`obm`,`bartlett`,`tukey`. `bm` represents batch
#'   means estimator, `obm` represents the overlapping batch means estimator,
#'   and `bartlett` and `tukey` represent the modified-Bartlett window and
#'   the Tukey-Hanning windows for the spectral variance estimators.
#' @param g A function that represents features of interest. g is applied to each row of x and
#'   thus g should take a vector input only. If g is NULL, g is set to be identity, which
#'   is estimation of the mean of the target density.
#'   
#' @return A value of the optimal batch size is returned.
#' 
#' @references 
#' Liu, Y., Vats, D., and Flegal, J. M. Batch size selection for variance estimators in MCMC, arXiv
#' preprint arXiv:1804.05975 (2019).
#' 
#' @seealso \code{\link{mcse.multi}}, which calls on batchSize. \code{\link{mcse}}, which calls on batchSize.
#' 
#' @export
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
#' batchSize(out)
#' batchSize(out, method = "obm")
#' batchSize(out, method = "bartlett")
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' X = multivariate_Gibbs_normal(n, mu, sigma)
#' batchSize(X)
#' batchSize(X, method = "obm")
#' batchSize(X, method = "bartlett")

batchSize <- function(x, method = c("bm", "obm", "bartlett", "tukey"), g = NULL)
{

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

  n <- dim(chain)[1]
  ar_fit <- apply(chain, 2, arp_approx)^2
  coeff <- ( sum(ar_fit[1,])/sum(ar_fit[2,]) )^(1/3)
  method = match.arg(method)

  b.const <- (3/2*n)*(method == "obm" || method == "bartlett" || method == "tukey") + (n)*(method == "bm")
  b <- b.const^(1/3) * coeff
  if(b <= 1) b <- 1

  b <- floor(b)
  return(b)
}


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

#' Multivariate Monte Carlo standard errors for expectations.
#' 
#' Function returns the estimate of the covariance matrix in the Markov Chain CLT using batch means
#' or spectral variance methods (with different lag windows). The function also returns the Monte
#' Carlo estimate.
#' 
#' @usage mcse.multi(x, method = "bm", r = 3, size = NULL, g = NULL, adjust = TRUE, blather = FALSE)
#' 
#' @param x A matrix or data frame of Markov chain output. Number of rows is the Monte
#' Carlo sample size.
#' @param method Any of `bm`,`obm`,`bartlett`,`tukey`. `bm` represents batch
#'   means estimator, `obm` represents the overlapping batch means estimator,
#'   and `bartlett` and `tukey` represent the modified-Bartlett window and
#'   the Tukey-Hanning windows for the spectral variance estimators.
#' @param r The lugsail parameter that converts a lag window into its lugsail equivalent.
#'   Larger values of `r` will typically imply less underestimation of ''cov'',
#'   but higher variability of the estimator. Default is `r = 3` and `r = 1,2` are
#'   good choices. `r > 5` is not recommended. Non-integer values are ok.
#' @param size Represents the batch size in "bm" and the truncation point in "bartlett" and
#'  "tukey". Default is NULL which implies that an optimal batch size is calculated
#'  using the batchSize() function. Can take character values of `sqroot` and
#'  `cuberoot` or any numeric value between 1 and n/2. `sqroot` means
#'  size is floor(n^(1/2)) and "cuberoot" means size is floor(n^(1/3)).
#' @param g A function that represents features of interest. g is applied to each row of x and
#'   thus g should take a vector input only. If g is NULL, g is set to be identity, which
#'   is estimation of the mean of the target density.
#' @param adjust Defaults to `TRUE`. logical for whether the matrix should automatically be ad-
#'  justed if unstable.
#' @param blather If `TRUE`, returns under-the-hood workings of the package.
#' 
#' @return A list is returned with the following components,
#' \describe{
#'  \item{`cov`}{a covariance matrix estimate.}
#'  \item{`est`}{estimate of g(x).}
#'  \item{`nsim`}{number of rows of the input x.}
#'  \item{`method`}{method used to calculate matrix cov.}
#'  \item{`size`}{value of size used to calculate cov.}
#'  \item{`adjust.used`}{whether an adjustment was used to calculate cov.}
#' }
#' @references 
#'  Vats, D., Flegal, J. M., and, Jones, G. L (2019) Multivariate Output Analysis for Markov chain
#'  Monte Carlo, Biometrika.
#'  Vats, D., Flegal, J. M., and, Jones, G. L. (2018) Strong Consistency of multivariate spectral   #'  variance estimators for Markov chain Monte Carlo, Bernoulli.
#'  
#' @seealso \code{\link{batchSize}}, which computes an optimal batch size. \code{\link{mcse.initseq}}, which computes
#' an initial sequence estimator. \code{\link{mcse}}, which acts on a vector. \code{\link{mcse.mat}}, which applies mcse
#' to each column of a matrix or data frame. \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard
#' errors for quantiles.
#'  
#' @export
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
#' mcse.bm <- mcse.multi(x = out)
#' mcse.tuk <- mcse.multi(x = out, method = "tukey")
#' # If we are only estimating the mean of the first component,
#' # and the second moment of the second component
#' g <- function(x) return(c(x[1], x[2]^2))
#' mcse <- mcse.multi(x = out, g = g)
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' X = multivariate_Gibbs_normal(n, mu, sigma)
#' mcse.bm <- mcse.multi(x = X)
#' mcse.tuk <- mcse.multi(x = X, method = "tukey")
#' # If we are only estimating the mean of the first component,
#' # and the second moment of the second component
#' g <- function(x) return(c(x[1], x[2]^2))
#' mcse <- mcse.multi(x = X, g = g)

mcse.multi <- function(x, method = c("bm", "obm", "bartlett", "tukey", "lug"), lug_params = c(3,0.5), size = NULL, g = NULL, adjust = TRUE, blather = FALSE)
{ 
  method = match.arg(method)
  
  # at some point the method used may be different
  # from the method asked. Blather will output this
  method.used <- method
  if(method == "lug")   # not releaved to the public. Inside option for developers
  {
    method = "bm"
    r <- 3
    c <- 0.5
  }
  
  r = lug_params[1]
  c = lug_params[2]
  
  if(!is.numeric(r)) stop("r should be numeric")
  if(!is.numeric(c)) stop("c should be numeric")
  
  if(r > 5) warning("We recommend using r <=5. r = 1,2,3 are standard")
  if(r < 1)  {
    warning("r cannot be less than 1. Setting r = 3")
    r = 3
  }
  if(c > 1) {
    warning("c cannot be greater than 1. Setting c = 0.5")
    c = 0.5
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
  n = dim(chain)[1]
  p = dim(chain)[2]

  # Setting batch size based on chosen option
  if(is.null(size))
  {
    b <- batchSize(x = x, method = method, g = g)  # optimal
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
      warning("size is either too large, too small, or not a number. Setting 'size' to the optimal value")
      size = batchSize(x = x, method = method, g = g)
    }
        

    b <- floor(size)
  }
  a <- floor(n/b)
  if(b == 1 && r != 1)
  {
    r = 1
  }
 ########################## 


  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)  # this is based on the full output. not n = a*b
  
  ## Setting matrix sizes to avoid dynamic memory 
  sig.mat = matrix(0, nrow = p, ncol = p)

  b = max(b, 2*r)

  message <- ""   # will store some info for blather
  
  ## Batch Means
  if(method == "bm")  {
    init.mat = mbmC(chain, b)
    sig.mat = init.mat
    if(r>1) {
      sig.mat <- (1/(1-c))*init.mat - (c/(1-c))*mbmC(chain, floor(b/r))
    }
  }
    
  
  ## Overlapping Batch Means
  if(method == "obm") {
    init.mat = mobmC(chain, b)
    sig.mat = init.mat
    if(r>1) {
      sig.mat <- (1/(1-c))*init.mat - (c/(1-c))*mobmC(chain, floor(b/r))
    }
  }
    
  
  ## Bartlett or Tukey
  if((method == "bartlett") || (method == "tukey"))
  {
    chain <- scale(chain, center = mu.hat, scale = FALSE)
    init.mat <-  mSVEfft(A = chain, b = b, method = method)
    sig.mat = init.mat
    if(r>1) {
      sig.mat <- (1/(1-c))*init.mat - (c/(1-c))*mSVEfft(A = chain, b = floor(b/r), method = method)
    }
  }
  
  method.used = paste("Lugsail ", method, " with r = ", r)
  if(prod(diag(sig.mat) > 0) == 0)  # If diagonals are negative, cannot use larger values of r
  {
    sig.mat <- init.mat
    method.used <- method
    message <- paste(message, paste("Diagonals were negative with r = ", r,". r = 1 was used.", sep = ""), sep = "")
  }
  
  sig.eigen <- eigen(sig.mat, only.values = TRUE)$values
  
  adjust.used <- FALSE  #whether an adjustment was made. None yet
  sig.adj = NULL
  
  if(adjust) # if adjust is FALSE, may output non PD estimator
  {
    if(min(sig.eigen) <= 0)  #needs an adjustment. No need to adjust is not needed
    {
      adjust.used <- TRUE
      warning("Estimated matrix not positive definite. The chain might be highly correlated or very high dimensional. Consider increasing the sample size. Using the default batch means estimator as a substitute.")
      if(method == "bm")  {
        sig.mat = init.mat
      } else  {
        sig.mat = mcse.multi(x, method = "bm", r = 1, size = size, g = g, adjust = FALSE, blather = FALSE)$cov
      }
      sig.adj = sig.mat
      sig.eigen <- eigen(sig.mat, only.values = TRUE)$values
    }    
  } 

  if(blather)
  {
    value = list("cov" = sig.mat, "est" = mu.hat, "nsim" = n, "eigen-values" = sig.eigen,
        "method" = method.used, "size" = b, "Adjustment-Used" = adjust.used, "message" = message,
        "cov.adj"=sig.adj)
  } else {
    value = list("cov" = sig.mat, "est" = mu.hat, "nsim" = n, "eigen-values" = sig.eigen,
                 "cov.adj"=sig.adj)
  }
  class(value) = "mcmcse"
  value

}


is.mcmcse <- function(x) inherits(x, "mcmcse")


