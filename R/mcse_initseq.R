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
#' @param adjust Logical; if \code{TRUE}, an adjustment is made to increase slightly the eigenvalues of
#'   the initial sequence estimator. The default is \code{FALSE}.
#' @param g A function that represents features of interest. \code{g} is applied to each row of x and
#'   thus \code{g} should take a vector input only. If \code{g} is \code{NULL}, \code{g} is set to be identity, which
#'   is estimation of the mean of the target density.
#' @param blather if \code{TRUE}, outputs under the hood information about the function.
#' 
#' @return A list is returned with the following components,
#'  \item{cov}{a covariance matrix estimate using intial sequence method.}
#'  \item{cov.adj}{a covariance matrix estimate using adjusted initial sequence method if the
#'  input \code{adjust=TRUE}.}
#'  \item{eigen_values}{eigen values of the estimate cov.}
#'   \item{method}{method used to calculate matrix cov.}
#'   \item{est}{estimate of g(x).}
#'   \item{nsim}{number of rows of the input x. Only if \code{blather = TRUE}.}
#'   \item{Adjustment_Used}{logical of whether an adjustment was made to the initial sequence estimator.
#'   Only if \code{blather = TRUE}.}
#' 
#' @references 
#' Dai, N and Jones, G.L. (2017)  Multivariate initial sequence estimators in Markov chain Monte Carlo, 
#' \emph{ Journal of Multivariate Analysis}, \bold{159}, 184-199.
#' 
#' @seealso 
#' \code{\link{mcse}}, which acts on a vector. 
#' \code{\link{mcse.mat}}, which applies \code{mcse} to each
#' column of a matrix or data frame. 
#' \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard errors for quantiles. 
#' \code{mcse.multi}, which estimates the covariance matrix in the
#' Markov Chain CLT using batch means or spectral variance methods.
#' 
#' @examples 
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu <- c(2, 50)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' out <- BVN_Gibbs(n, mu, sigma)
#' 
#' out.mcse <- mcse.initseq(x = out)
#' out.mcse.adj <- mcse.initseq(x = out, adjust = TRUE)
#' 
#' # If we are only estimating the mean of the first component,
#' # and the second moment of the second component
#' g <- function(x) return(c(x[1], x[2]^2))
#' out.g.mcse <- mcse.initseq(x = out, g = g)
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
  
  if(n < (p+1))
    stop("sample size is insufficient for a Markov chain of this dimension")
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
    value = list("cov" = sig, "cov.adj"=sig.adj, "method" = "initial sequence", "eigen_values" =                    sig.eigen, "est" = mu.hat, "nsim" = n, "Adjustment_Used" = adjust)
  } else{
    value = list("cov" = sig, "cov.adj"=sig.adj, "nsim" = n, "eigen_values" = sig.eigen,
                 "est" = mu.hat)
  }
  class(value) = "mcmcse"
  value
}



