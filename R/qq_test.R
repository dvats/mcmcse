# qqTest <- function(x, covmat, g = NULL)
# {
# 	chain <- as.matrix(x)
# 	if(!is.matrix(x) && !is.data.frame(x))
# 	  stop("'x' must be a matrix or data frame.")

# 	if (is.function(g)) 
# 	{
# 	  chain <- apply(x, 1, g)

# 	  if(is.vector(chain))
# 	  {
# 	    chain <- as.matrix(chain)
# 	  }else
# 	  {
# 	    chain <- t(chain)
# 	  }
# 	}

# 	mu <- colMeans(chain)
# 	n <- dim(chain)[1]
# 	p <- dim(chain)[2]
# 	decomp  <- svd(covmat)
# 	inv.root <- decomp$v %*% diag( (decomp$d^(-1/2)), p) %*% t(decomp$u)
# 	qqnorm(inv.root%*%mu)
# 	qqline(inv.root%*%mu)
# }

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