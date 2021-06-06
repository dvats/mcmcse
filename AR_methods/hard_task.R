##################################################
## Aim: Implement effocoent RWMG for mvt Gaussian.
## Also: Calculate numerically stable ESS. 
##################################################
set.seed(10)

# Mean is 0 and covariance matrix is identity
log_unnormalised_posterior <- function(X)  {
  ans = (-0.5) * sum(X^2)
}

# Random walk MH
RWMH <- function(n, init, h)  {
    p = length(init)
    accept = 0
    output = matrix(, nrow = n, ncol = p)
    output[1,] = init
    h_root = sqrt(h)
    for(t in 2:n) {
      prop =  output[t-1,] + rnorm(p,0,h_root) 
      log_ratio = log_unnormalised_posterior(prop) - log_unnormalised_posterior(output[t-1,]) # work with log for numerical stability
      if(log(runif(1)) < log_ratio) {
        output[t,] = prop
        accept = accept + 1
      }
      else  {
        output[t,] = output[t-1,]
      }
    }
    print(accept/n) # acceptance probability
    output
}

# Divide X into a batches of b length each
make_batch <- function(X, n, a, b)  {
  batches = matrix(, nrow = a, ncol = dim(X)[2])
  
  for(t in a:1) # removing the initial samples rather than the final samples if n != ab
    batches[t,] = colMeans(X[(b*(t-1)+1):(b*t),])
  
  batches
}

# Estimate CLT covariance matrix
bm_estimator <- function(chain) {
  n = dim(chain)[1]
  p = dim(chain)[2]
  a = floor(sqrt(n)) # no of batches is sqroot(n)
  b = a
  X = make_batch(chain, n, a, b) # divide the chain into batches
  chain_mean = colMeans(X)
  bm = (t(X - chain_mean) %*% (X - chain_mean)) * (b/(a-1)) # formula for batch means from Vats et al 2017
  bm
}

# Modify eigen values if they are below epsilon
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

# Estimate effective sample size
ess <- function(chain)  {
  n = dim(chain)[1]
  p = dim(chain)[2]
  sample_var = cov(chain) # Sample covariance matrix
  CLT_var = bm_estimator(chain) # CLT covariance estimate
  if (min(eigen(CLT_var, only.values = TRUE)$values) <= 0)  {
    CLT_var = adjust_matrix(CLT_var, N=n) # Tweak the matrix if BM estimate is not pd
  }
  sample_var_det = exp(sum(log(eigen(sample_var, only.values = TRUE)$values))/p) # Det ^ 1/p
  CLT_var_det = exp(sum(log(eigen(CLT_var, only.values = TRUE)$values))/p) # Using exp(log(sum(eigen))) for numerical stability.
  ess = n*(sample_var_det/CLT_var_det) # formula for ESS from Vats et al 2017
  ess
}

n = 1e4
p = 100
init = numeric(p)
h = 0.05
ans = RWMH(n, init, h)
ess(ans)
