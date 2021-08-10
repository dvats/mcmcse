library(testthat)
library(mcmcse)

context("validResults")

test_mcse_multi <- function(x, method = "bm", size = "sqroot", g = NULL, 
                     level = 0.95, large = FALSE)
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
  
  ## Setting dimensions on the mcmc output. 
  n = dim(chain)[1]
  p = dim(chain)[2]
  
  
  ## Initializing b_n 
  if(size == "sqroot")
  {
    b = floor(sqrt(n))
  } 
  else if(size == "cuberoot") {
    b = floor(n^(1/3))
  }
  else {
    if (!is.numeric(size) || size <= 1 || size > n) 
      stop("'size' must be a numeric quantity not larger than n.")
    b = floor(size)
  }
  a = floor(n/b)
  m <- 0
  
  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)
  
  ## Batch means
  batch_means <- matrix(0, nrow = a, ncol = p)
  
  ## Setting matrix sizes to avoid dynamic memory 
  sig.mat = matrix(0, nrow = p, ncol = p)
  
  ## only does bm and obm
  if(method != "bm" && method != "bartlett" && method != "tukey")
  {
    stop("No such method available")
  }
  
  ## Batch Means
  if(method == "bm")
  {
    
    for(i in 1:a)
    {
      batch_means[i,] = colMeans(chain[((i-1)*b+1): (i*b),])
    }
    
    for(j in 1:a)
    {
      sig.mat <- sig.mat + tcrossprod(batch_means[j,] - mu.hat)
    }
    
    sig.mat <- b*sig.mat/(a-1)
    m <- a - 1
  }
  
  ## Modified Bartlett Window
  
  if(method == "bartlett" || method == "tukey")
  {
    m <- n - b
    lag <- function(x)
    { 
      ifelse(method == "bartlett",  (1 - x/b), (1 + cos(pi * x/b))/2 ) 
    }
    
    ## centering
    for(i in 1:n)
    {
      chain[i,] <- chain[i,] - mu.hat
    }
    
    for(s in 1:(b-1))
    {
      for(j in 1: (n-s))
      {
        foo <- chain[j,]%*%t(chain[j+s,])
        sig.mat <- sig.mat + lag(s)*(foo + t(foo))
      }
    }
    sig.mat <- (sig.mat + t(chain)%*%chain)/n
  }
  
  crit <- ifelse(large, qchisq(level, df = p)/n,
                 exp(log(p) + log(m) - log(n) - log(m-p+1) + log(qf(level, p, m-p+1))) )
  dummy <- log(2) + (p/2)*log(pi*crit) - log(p) - lgamma(p/2)
  # log.det.sig <- sum(log(eigen(sig.mat, only.values = TRUE)$values))
  det.sig <- det(sig.mat)
  volume.sig <- exp(dummy + (1/2)*log(det.sig))
  pth.vol <- (volume.sig)^(1/p)
  return(list("cov" = sig.mat, "vol" = pth.vol, "est" = mu.hat, "nsim" = n, 
              "method" = method, "large" = large, "size" = b))
}

n <- 1e3
p <- 3
b <- floor(sqrt(n))
out <- matrix(rnorm(n*p), nrow = n, ncol = p)

mbatch <- test_mcse_multi(out)
mbart <- test_mcse_multi(out, method = "bartlett")
mtukey <- test_mcse_multi(out, method = "tukey")

test_that("test if functions return correct values for mcse.multi", {
  
  expect_equal(mcse.multi(out), mbatch)
  expect_equal(mcse.multi(out, method = "bartlett"), mbart)
  expect_equal(mcse.multi(out, method = "tukey"), mtukey)
  
})

test_mcse_q = function(x, q, size = "sqroot", g = NULL, method = c("bm", "obm", "sub"), warn = FALSE)
{
  if (! is.function(g))
    g = function(x) return(x)
  n = length(x)
  if (n < 1000)
  {
    if (warn)
      warning("too few samples (less than 1,000)")
    if (n < 10)
      return(NA)
  }
  if (size == "sqroot") 
  {
    b = floor(sqrt(n))
    a = floor(n / b)
  }
  else
  {
    if (! is.numeric(size) || size < 1 || size == Inf)  {
      warning("'size' must be a finite numeric quantity larger than 1. Setting 'size' to the optimal value")
      size = batchSize(x = x, method = method, g = g)
    }
    
    b = floor(size)
    a = floor(n / b)
  }
  method = match.arg(method)
  counting = function(var.vector, var.number)
  {
    return(length(var.vector[var.vector <= var.number]))
  }
  if (! is.numeric(q) || q <= 0 || q >= 1)    {
    stop("'q' must be from (0, 1).")
  }
  
  quant = function(input) { quantile(input, prob = q, type = 1, names = FALSE) }
  if (method == "bm")
  {
    xi.hat = quant(g(x))
    y = sapply(1:a, function(k) return(counting(g(x[((k - 1) * b + 1):(k * b)]), xi.hat))) / b
    mu.hat = mean(y)
    var.hat = b * sum((y - mu.hat)^2) / (a - 1)
    f.hat.junk = density(g(x), from = xi.hat, to = xi.hat, n = 1)
    f.hat = f.hat.junk$y
    se = sqrt(var.hat / n) / f.hat
  }
  else if (method == "obm")
  {
    xi.hat = quant(g(x))
    a = n - b + 1
    y = sapply(1:a, function(k) return(counting(g(x[k:(k + b - 1)]), xi.hat))) / b
    mu.hat = mean(y)
    var.hat = n * b * sum((y - mu.hat)^2) / (a - 1) / a
    f.hat.junk = density(g(x), from = xi.hat, to = xi.hat, n = 1)
    f.hat = f.hat.junk$y
    se = sqrt(var.hat / n) / f.hat
  } 
  else # method == "sub"
  {
    xi.hat = quant(g(x))
    a = n - b + 1
    y = sapply(1:a, function(k) return(quant(g(x[k:(k + b - 1)]))))
    mu.hat = mean(y)
    var.hat = n * b * sum((y - mu.hat)^2) / (a - 1) / a
    se = sqrt(var.hat / n)
  }
  value = list("est" = xi.hat, "se" = se, "nsim" = n)
  class(value) = "mcmcse"
  value
}

n = 1e4
out = double(n)
out[1] = 2
for (i in 1:(n - 1))
  out[i + 1] = 0.9 * out[i] + rnorm(1)

mbatch <- test_mcse_q(out)
mobm <- test_mcse_q(out, method = "obm")
msub <- test_mcse_q(out, method = "sub")

test_that("test if functions return correct values for mcse.multi", {
  
  expect_equal(mcse.q(out), mbatch)
  expect_equal(mcse.q(out, method = "obm"), mobm)
  expect_equal(mcse.q(out, method = "sub"), msub)
  
})


Gibbs_sampler <- function(mu1, mu2, a, b, rho, init, n) {
  X <- matrix(0, nrow = n, ncol = 2)
  X[1, 1] = init[1]
  X[1, 2] = init[2] 
  for (i in 2:n) {
    X[i, 1] = rnorm(1, mu1 + (rho / b) * (X[i - 1, 2] - mu2), sqrt(a - (rho ^ 2) / b))
    X[i, 2] = rnorm(1, mu2 + (rho / a) * (X[i, 1] - mu1), sqrt(b - (rho ^ 2) / a))
  }
  return(X)
}


test_that("ess performs univariate sampling and output makes sense", {
  chain = Gibbs_sampler(2, 50, 1, 1, 0.5, c(2, 50), 1e4)
  ess_est = ess(chain)
  expect_equal(length(ess_est), 2)
  sapply(1:2, function(i) expect_lte(ess_est[i], 1e4)) # not right
})
