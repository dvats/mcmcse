multivariate_Gibbs_normal <- function(n, mu, sigma) {
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

