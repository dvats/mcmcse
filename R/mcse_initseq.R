# sourceCpp("inseq.cpp")
#' Multivariate Monte Carlo standard errors for expectations with the initial sequence method of Dai
#' and Jones (2017)
#' 
#' 
#' Function returns the estimate of the covariance matrix in the Markov Chain central limit theorem
#' using initial sequence method. This method is designed to give an asymptotically conservative
#' estimate of the Monte Carlo standard error.
#' 
#' @usage mcse.initseq(x, g = NULL, adjust = FALSE, blather = FALSE)
#' 
#' @param x A matrix or data frame of Markov chain output. Number of rows is the Monte
#'   Carlo sample size.
#' @param adjust Logical; if `TRUE`, an adjustment is made to increase slightly the eigenvalues of
#'   the initial sequence estimator. The default is `FALSE`.
#' @param g A function that represents features of interest. g is applied to each row of x and
#'   thus g should take a vector input only. If g is `NULL`, g is set to be identity, which
#'   is estimation of the mean of the target density.
#' @param blather if `TRUE`, outputs under the hood information about the function.
#' 
#' @return A list is returned with the following components,
#' \describe{
#'  \item{`cov`}{a covariance matrix estimate using intial sequence method}
#'  \item{`cov.adj`}{a covariance matrix estimate using adjusted initial sequence method if the
#'  input `adjust=TRUE.`}
#'   \item{`est`}{estimate of g(x).}
#'   \item{`nsim`}{number of rows of the input x. Only if `blather = TRUE`.}
#'   \item{`adjust`}{logical of whether an adjustment was made to the initial sequence estimator.
#'   Only if `blather = TRUE`.}
#' }
#' 
#' @references 
#' Dai, N and Jones, G.L. (2017+) Multivariate initial sequence estimators in Markov chain Monte
#' Carlo, Journal of Multivariate Analysis.
#' 
#' @seealso `initseq{mcmc}`, which is a different univariate initial sequence estimator.
#' \code{\link{mcse}}, which acts on a vector. \code{\link{mcse.mat}}, which applies `mcse` to each
#' column of a matrix or data frame. \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which
#' compute standard errors for quantiles. mcse.multi, which estimates the covariance matrix in the
#' Markov Chain CLT using batch means or spectral variance methods.
#' 
#' @examples 
#' library(mAr)
#' p <- 3
#' n <- 1000
#' omega <- 5*diag(1,p)
#' ## Making correlation matrix var(1) model
#' set.seed(100)
#' foo <- matrix(rnorm(p^2), nrow = p)
#' foo <- foo %*% t(foo)
#' phi <- foo / (max(eigen(foo)$values) + 1)
#' dat <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))
#' out.mcse <- mcse.initseq(x = dat)
#' out.mcse.adj <- mcse.initseq(x = dat,adjust = TRUE)
#' # If we are only estimating the mean of the first component,
#' # and the second moment of the second component
#' g <- function(x) return(c(x[1], x[2]^2))
#' out.g.mcse <- mcse.initseq(x = dat, g = g)
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' X = multivariate_Gibbs_normal(n, mu, sigma)
#' out.mcse <- mcse.initseq(x = X)
#' out.mcse.adj <- mcse.initseq(x = X,adjust = TRUE)
#' # If we are only estimating the mean of the first component,
#' # and the second moment of the second component
#' g <- function(x) return(c(x[1], x[2]^2))
#' out.g.mcse <- mcse.initseq(x = X, g = g)
#' 
#' @export
#' 

mcse.initseq <- function(x, g = NULL, adjust = FALSE, blather = FALSE)
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
  
  ## Setting dimensions on the mcmc output
  n = dim(chain)[1]
  p = dim(chain)[2]
  
  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)

  ## Initial Sequence Estimator(s)
  res <- inseq(chain, adjust)

  ##sig=initial sequence estimator without asjustment
  sig <- res$Sig
  
  ##sig.adj=initial sequence estimator with asjustment, if adjust=T
  ##       =NULL, if adjust=F
  if(adjust)
  {
    sig.adj <- res$Sigadj
  }else
  {
    sig.adj <- NULL
  }

  sig.eigen = eigen(sig, only.values = TRUE)$values
  if(blather)
  {
    value = list("cov" = sig, "cov.adj"=sig.adj, "method" = "initial sequence", "eigen-values" =                    sig.eigen, "est" = mu.hat, "nsim" = n, "Adjustment-Used" = adjust, "size" = NULL,
                 message = "")
  } else{
    value = list("cov" = sig, "cov.adj"=sig.adj, "nsim" = n, "eigen-values" = sig.eigen,
              "est" = mu.hat)
  }
  class(value) = "mcmcse"
  value
}



