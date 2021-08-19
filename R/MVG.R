#' MCMC on a bivariate normal distribution 
#' 
#' Function returns Gibbs samples from a bivariate normal target density.
#' 
#' 
#' @usage BVN_Gibbs(n, mu, sigma)
#' 
#' @param n Sample size of the Markov chain.
#' @param mu A 2 dimensional vector. Mean of the target normal distribution.
#' @param sigma 2 x 2 symmetric positive semi-definite matrix. The covariance matrix of the target normal distribution.
#'   
#' @return An n x 2 matrix of the Gibbs samples.
#' 
#' @export
#' 
#' @examples 
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' X = BVN_Gibbs(n, mu, sigma)
#' 

BVN_Gibbs <- function(n, mu, sigma) {
  init <- matrix(rnorm(ncol(sigma)), nrow = 1, byrow = TRUE) %*% chol(sigma) + mu
  rho = sigma[1,2]
  a = sigma[1,1]
  b = sigma[2,2]
  X <- matrix(0, nrow = n, ncol = 2)
  X[1, ] = init
  ## Gibbs sampler to generate the Markov chain
  for (i in 2:n) {
    X[i, 1] = rnorm(1, mu[1] + (rho / b) * (X[i - 1, 2] - mu[2]), sqrt(a - (rho ^ 2) / b))
    X[i, 2] = rnorm(1, mu[2] + (rho / a) * (X[i, 1] - mu[1]), sqrt(b - (rho ^ 2) / a))
  }
  X
}

