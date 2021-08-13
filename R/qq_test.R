
#' QQplot for Markov chains
#' 
#' QQplot for Markov chains using an estimate of the Markov Chain CLT covariance matrix.
#' 
#' @usage qqTest(mcse.obj)
#' 
#' @param mcse.obj the list returned by the `mcse.multi` or `mcse.initseq` command
#' 
#' @examples 
#' library(mAr)
#' p <- 35
#' n <- 1e4
#' omega <- 5*diag(1,p)
#' ## Making correlation matrix var(1) model
#' set.seed(100)
#' foo <- matrix(rnorm(p^2), nrow = p)
#' foo <- foo %*% t(foo)
#' phi <- foo / (max(eigen(foo)$values) + 1)
#' out <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))
#' mcse.bm <- mcse.multi(x = out)
#' qqTest(mcse.bm)
#' mcse.isadj <- mcse.initseq(x = out, adjust = TRUE)
#' qqTest(mcse.isadj)
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' X = BVN_Gibbs(n, mu, sigma)
#' mcse.bm <- mcse.multi(x = X)
#' qqTest(mcse.bm)
#' mcse.isadj <- mcse.initseq(x = X, adjust = TRUE)
#' qqTest(mcse.isadj)
#' 
#' @export
#' 

qqTest <- function(mcse.obj)
{
  mu <- mcse.obj$est
  n <- mcse.obj$nsim
  p <- length(mcse.obj$est)
  
  if(sum(names(mcse.obj) == "adjust"))
  {
    if(mcse.obj$adjust)
    {
      covmat <- mcse.obj$cov.adj
    }else{
      covmat <- mcse.obj$cov
    }
  }else{
    covmat <- mcse.obj$cov
  }
  decomp  <- svd(covmat)
  inv.root <- decomp$v %*% diag( (decomp$d^(-1/2)), p) %*% t(decomp$u)
  qqnorm(inv.root%*%mu)
  qqline(inv.root%*%mu)
}
