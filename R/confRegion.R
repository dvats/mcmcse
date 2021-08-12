library(ellipse)
#' Confidence regions (ellipses) for Monte Carlo estimates
#' 
#' Constructs confidence regions (ellipses) from the Markov chain output for the features of interest
#' Function uses the ellipse package.
#' 
#' @usage confRegion(mcse.obj, which = c(1,2), level = .95)
#' 
#' @param mcse.obj The list returned by the mcse.multi or mcse.initseq command.
#' @param which Integer vector of length 2 indicating the component for which to make the confidence
#' ellipse. Chooses the first two by #'  default.
#' @param level confidence level for the ellipse.
#' 
#' @return Returns a matrix of x and y coordinates for the ellipse. Use plot function on the matrix
#' to plot the ellipse.
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
#' mcerror <- mcse.multi(out, blather = TRUE)
#' ## Plotting the ellipse
#' plot(confRegion(mcerror), type = 'l')
#' 
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' X = multivariate_Gibbs_normal(n, mu, sigma)
#' mcerror <- mcse.multi(X, blather = TRUE)
#' ## Plotting the ellipse
#' plot(confRegion(mcerror), type = 'l')
#' 

confRegion <- function(mcse.obj, which = c(1,2), level = .95)
{
  mat <- mcse.obj$cov
  n <- mcse.obj$nsim
  p <- 2
  
  crit <- qchisq(level, df = p)/n
  
  mu <- mcse.obj$est
  
  return(ellipse(mat, centre = mu[which], t = sqrt(crit), which = which))       
  
}