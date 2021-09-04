set.seed(10)
library(testthat)
library(mcmcse)


context("validResults")

eureka <- function(order_max, r, g, coefs, var, a, threshold) {
  v = r[1]
  d = r[2]
  a[1] = 1.0
  coefs[1,1] = g[2]/v
  
  ans = list("vars" = var, "coefs" = coefs, "order" = order_max)
  
  if(abs(coefs[1,1]) <= threshold) {
    ans$vars = var
    ans$coefs = coefs
    ans$order = 0
    return(ans)
  }
  
  q = coefs[1,1] * r[2]
  var[1] = (1 - coefs[1,1]*coefs[1,1]) * r[1]
  
  if(order_max == 1) {
    ans$vars = var
    ans$coefs = coefs
    ans$order = 1
    return(ans)
  }
  
  for(l in 2:order_max) {
    a[l] = -d/v
    
    if(l > 2) {
      l1 = (l-2)/2
      l2 = l1+1
      
      for(j in 2:l2) {
        hold = a[j]
        k = l-j+1
        a[j] = a[j] + a[l] * a[k]
        a[k] = a[k] + a[l] * hold
      }
      
      if((2*l1) != (l-2))
        a[l2+1] = a[l2+1] * (1.0 + a[l])
    }
    
    v = v + a[l] * d
    coefs[l,l] = (g[l+1] - q)/v
    
    if(abs(coefs[l,l]) <= threshold) {
      ans$vars = var
      ans$coefs = coefs
      ans$order = l-1
      return(ans)
    }
    
    
    for(j in 1:(l-1)) {
      coefs[l,j] = coefs[l-1,j] + coefs[l,l] * a[l-j+1]
    }
    
    var[l] = var[l-1] * (1 - coefs[l,l]*coefs[l,l])
    
    if(l == order_max) {
      ans$vars = var
      ans$vars = coefs
      ans$order = order_max
      return(ans)
    }
    
    d = 0.0
    q = 0.0
    
    for(i in 1:l) {
      k = l-i+2
      d = d + a[i] * r[k]
      q = q + coefs[l,i] * r[k]
    }
    
  }
  ans$vars = var
  ans$coefs = coefs
  ans$order = order_max
  return(ans)
}

ar_yw <- function(order.max, r, g, coefs, var, a, threshold, n) {
  ans = eureka(order.max, r, g, coefs, var, a, threshold)
  coef_vec = numeric(ans$order)
  
  if(ans$order> 0)  {
    coef_vec = t(ans$coefs[ans$order, 1:ans$order])
    var_pred = ans$vars[ans$order]
  }
  else{
    var_pred = r[1]
  }
  
  var_pred = var_pred*n/(n-(ans$order+1))
  
  ret = list("coefs" = coef_vec, "vars" = var_pred, "order" = ans$order)
  return(ret)
}

arp_approx <- function(xacf, order.max, n)
{
  threshold = qnorm((1.95)/2)/sqrt(n)
  # Fitting a univariate AR(m) model
  ar.fit <- ar_yw(order.max, xacf, xacf, matrix(0, nrow = order.max, ncol = order.max), 
                  numeric(order.max), numeric(order.max), threshold, n)
  
  # estimated autocovariances
  gammas <- xacf #as.numeric(acf(x, type = "covariance", lag.max = ar.fit$order, plot = FALSE)$acf)
  spec <- ar.fit$vars/(1-sum(ar.fit$coefs))^2  #asym variance
  
  if(ar.fit$order != 0)
  {
    foo <- 0
    for(i in 1:ar.fit$order)
    {
      for(k in 1:i)
      {
        foo <- foo + ar.fit$coefs[i]*k*gammas[abs(k-i)+1]
      }
    }
    Gamma <- 2*(foo + (spec - gammas[1])/2 *sum(1:ar.fit$order * ar.fit$coefs)  )/(1-sum(ar.fit$coefs))
  } else{
    Gamma <- 0
  }
  rtn <- cbind(Gamma, spec)
  colnames(rtn) <- c("Gamma", "Sigma")
  return(rtn)
}

test_batchSize <- function(x, method = c("bm", "obm", "bartlett", "tukey"), g = NULL)  {
  method = match.arg(method)
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
  p <- dim(chain)[2]
  
  if(n < (p+1))
    stop("sample size is insufficient for a Markov chain of this dimension")
  
  order.max <- min(p, n - 1L, floor(10 * log10(n))) # Maximum order up to which AR is fit
  xacf = matrix(, nrow = order.max+1, ncol = p)
  xacf = sapply(1:p, function(i) acf(chain[,i], type = "covariance", lag.max = order.max, plot =
                                       FALSE, demean=TRUE, na.action = na.pass)$acf)
  
  ar_fit <- apply(xacf, 2, arp_approx, order.max = order.max, n = n)^2
  coeff <- ( sum(ar_fit[1,])/sum(ar_fit[2,]) )^(1/3)
  method = match.arg(method)
  
  b.const <- (3/2*n)*(method == "obm" || method == "bartlett" || method == "tukey") + (n)*(method == "bm")
  b <- b.const^(1/3) * coeff
  if(b <= 1) b <- 1
  
  b <- floor(b)
  b <- min(b, floor(n / (p + 1)))
  if(n > 10)
    b = min(b, floor(n/10))
  b <- floor(b)
  return(b)
}

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

test_that("batchsize test", {
  mvg = Gibbs_sampler(2, 50, 1, 1, 0.5, c(2, 50), 1e4)
  expect_equal(test_batchSize(mvg), batchSize(mvg, fast = FALSE))
  expect_equal(test_batchSize(mvg, method = "obm"), batchSize(mvg, method = "obm", fast = FALSE))
  expect_equal(test_batchSize(mvg, method = "bartlett"), batchSize(mvg, method = "bartlett", 
                                                                   fast = FALSE))
  expect_equal(test_batchSize(mvg, method = "tukey"), batchSize(mvg, method = "tukey", 
                                                                fast = FALSE))
})

test_mcse_multi <- function(x, method = c("bm", "obm", "bartlett", "tukey", "lug"), r=3, size = NULL, g = NULL, adjust = TRUE, blather = FALSE)
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
  r=1
  c = 0.5
  
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
  
  if(n < (p+1))
    stop("sample size is insufficient for a Markov chain of this dimension")
  
  # Setting batch size based on chosen option
  if(is.null(size))
  {
    b <- batchSize(x = x, method = method)  # optimal
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
      warning("size is either too large, too small, or not a number. Setting 'size' to the optimal
              value")
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
  
  message <- ""   # will store some info for blather
  
  if(b < (2*r)) {
    r = 1       
    message = paste(message, "estimated batch size is low, lugsail not required")
  }

  ## Setting matrix sizes to avoid dynamic memory 
  sig.mat = matrix(0, nrow = p, ncol = p)
  sig.sum = matrix(0, nrow = p, ncol = p)
  
  ## only does bm and obm
  if(method != "bm" && method != "obm" && method != "bartlett" && method != "tukey")
  {
    stop("No such method available")
  }
  
  if(b == 1)
    sig.mat = var(chain)
  else  {
    if(method == "bm")
    {
      a = floor(n/b)
      y.bar <- matrix(0,nrow = a, ncol = p)
      y.bar <- apply(chain, 2, function(x) sapply(1:a, function(k) mean(x[((k-1)*b+1):(k*b)])))
      for(i in 1:a)
      {
        sig.sum = sig.sum + tcrossprod(y.bar[i,] - mu.hat)
      }
      
      sig.mat <- b*sig.sum/(a-1)
    }
    
    
    if(method == "obm")
    {
      a = n-b+1
      y.bar <- matrix(0,nrow = a, ncol = p)
      y.bar <- apply(chain, 2, function(x) sapply(1:a, function(k) mean(x[((k-1)*b+1):(k*b)])))
      y.bar <- apply(chain, 2, function(x) sapply(1:a, function(k) mean(x[k:(k+b-1)])))
      for(i in 1:a)
      {
        sig.sum = sig.sum + tcrossprod(y.bar[i,] - mu.hat)
      }
      
      sig.mat <- b*sig.sum/n
      
    }
    
    ## Modified Bartlett Window
    
    if(method == "bartlett")
    {
      chain <- scale(chain, center = mu.hat, scale = FALSE)
      dummy <- lapply(1:(b-1), function(j) (1 - j/b)*(  t(chain[1:(n-j), ])%*%chain[(j+1):n, ] 
                                                        + t(chain[(1+j):(n), ])%*%chain[1:(n-j), ]  ) )
      
      sig.sum <- t(chain)%*%chain + Reduce('+', dummy)
      
      sig.mat <- sig.sum/n
    }
    
    if(method == "tukey")
    {
      chain <- scale(chain, center = mu.hat, scale = FALSE)
      dummy <- lapply(1:(b-1), function(j) ((.5)*(1 + cos(pi * j/b))*t(chain[1:(n-j), ]))%*%chain[(j+1):n, ] 
                      + ((.5)*(1 + cos(pi * j/b))*t(chain[(1+j):(n), ]))%*%chain[1:(n-j), ])
      
      sig.sum <- t(chain[1:(n), ])%*%chain[(1):n, ] + Reduce('+', dummy)
      sig.mat <- sig.sum/n
      
    }
  }
  
  
  
  # ## Batch Means
  # if(method == "bm")
  # {
  # 
  #   for(i in 1:a)
  #   {
  #     batch_means[i,] = colMeans(chain[((i-1)*b+1): (i*b),])
  #   }
  # 
  #   for(j in 1:a)
  #   {
  #     sig.mat <- sig.mat + tcrossprod(batch_means[j,] - mu.hat)
  #   }
  # 
  #   sig.mat <- b*sig.mat/(a-1)
  #   m <- a - 1
  # }

  # ## Modified Bartlett Window
  # 
  # if(method == "bartlett" || method == "tukey")
  # {
  #   m <- n - b
  #   lag <- function(x)
  #   { 
  #     ifelse(method == "bartlett",  (1 - x/b), (1 + cos(pi * x/b))/2 ) 
  #   }
  #   
  #   ## centering
  #   for(i in 1:n)
  #   {
  #     chain[i,] <- chain[i,] - mu.hat
  #   }
  #   
  #   for(s in 1:(b-1))
  #   {
  #     for(j in 1: (n-s))
  #     {
  #       foo <- chain[j,]%*%t(chain[j+s,])
  #       sig.mat <- sig.mat + lag(s)*(foo + t(foo))
  #     }
  #   }
  #   sig.mat <- (sig.mat + t(chain)%*%chain)/n
  # }
  
  sig.eigen = eigen(sig.mat, only.values = TRUE)$values

  value = list("cov" = sig.mat, "est" = mu.hat, "nsim" = n, "eigen_values" = sig.eigen)
  class(value) = "mcmcse"
  value
}

n <- 1e3
p <- 3
b <- floor(sqrt(n))
out <- matrix(rnorm(n*p), nrow = n, ncol = p)

mbatch <- test_mcse_multi(out)
mobm <- test_mcse_multi(out, method = "obm")
mbart <- test_mcse_multi(out, method = "bartlett")
mtukey <- test_mcse_multi(out, method = "tukey")

test_that("test if functions return correct values for mcse.multi", {
  
  expect_equal(mcse.multi(out, r=1, adjust = FALSE), mbatch)
  expect_equal(mcse.multi(out, method = "obm", r=1, adjust = FALSE), mobm)
  expect_equal(mcse.multi(out, method = "bartlett", r=1, adjust = FALSE), mbart)
  expect_equal(mcse.multi(out, method = "tukey", r=1, adjust = FALSE), mtukey)
  
})


test_mcse.q = function(x, q, size = NULL, g = NULL, method = c("bm", "obm", "sub"), warn = FALSE)
{
  method = match.arg(method)
  if(is.function(g))
    x = sapply(x, g)
  
  n = length(x)
  
  if (n < 1000)
  {
    if (warn)
      warning("too few samples (less than 1,000)")
    if (n < 10)
      return(NA)
  }
  
  if(is.null(size)) {
    b <- batchSize(x = x, method = method)
  } else if(size == "sqroot")  {
    b <- floor(sqrt(n))
  } else if(size == "cuberoot") {
    b <- floor(n^(1/3))
  } else  {
    if (!is.numeric(size) || size < 1 || size >= n || floor(n/size) <=1) {
      warning("size is either too large, too small, or not a number. Setting 'size' to n^(1/2)")
      size = sqrt(n)
    }
    
    b <- floor(size)
  }
  
  a <- floor(n/b)
  
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
    xi.hat = quant(x)
    y = sapply(1:a, function(k) return(counting(x[((k - 1) * b + 1):(k * b)], xi.hat))) / b
    mu.hat = mean(y)
    var.hat = b * sum((y - mu.hat)^2) / (a - 1)
    f.hat.junk = density(x, from = xi.hat, to = xi.hat, n = 1)
    f.hat = f.hat.junk$y
    se = sqrt(var.hat / n) / f.hat
  }
  else if (method == "obm")
  {
    xi.hat = quant(x)
    a = n - b + 1
    y = sapply(1:a, function(k) return(counting(x[k:(k + b - 1)], xi.hat))) / b
    mu.hat = mean(y)
    var.hat = n * b * sum((y - mu.hat)^2) / (a - 1) / a
    f.hat.junk = density(x, from = xi.hat, to = xi.hat, n = 1)
    f.hat = f.hat.junk$y
    se = sqrt(var.hat / n) / f.hat
  } 
  else # method == "sub"
  {
    xi.hat = quant(x)
    a = n - b + 1
    y = sapply(1:a, function(k) return(quant(x[k:(k + b - 1)])))
    mu.hat = mean(y)
    var.hat = n * b * sum((y - mu.hat)^2) / (a - 1) / a
    se = sqrt(var.hat / n)
  }
  value = list("est" = xi.hat, "se" = se, "nsim" = n)
  value
}

n = 1e4
out = double(n)
out[1] = 2
for (i in 1:(n - 1))
  out[i + 1] = 0.9 * out[i] + rnorm(1)

mbatch <- test_mcse.q(out, q = 0.7)
mobm <- test_mcse.q(out, q = 0.7, method = "obm")
msub <- test_mcse.q(out, q = 0.7, method = "sub")

test_that("test if functions return correct values for mcse.multi", {
  
  expect_equal(mcse.q(out, q = 0.7), mbatch)
  expect_equal(mcse.q(out, q = 0.7, method = "obm"), mobm)
  expect_equal(mcse.q(out, q = 0.7, method = "sub"), msub)
  
})





test_that("ess performs univariate sampling and output makes sense", {
  chain = Gibbs_sampler(2, 50, 1, 1, 0.5, c(2, 50), 1e4)
  ess_est = ess(chain)
  expect_equal(length(ess_est), 2)
})

test_that("mcse performs univariate output analysis",{
  chain = Gibbs_sampler(2, 50, 1, 1, 0.5, c(2, 50), 1e4)
  foo = mcse.mat(chain)
  expect_equal(length(foo[,1]), 2)
})

# test_that("batchsize follows constraints", {
#   chain = matrix(rnorm(4), nrow = 2)
#   expect_error(batchSize(chain), "sample size is insufficient for a Markov chain of this dimension")
#   
#   x = numeric(20)
#   x[1] = 1
#   for(i in 2:20)
#     x[i] = 0.99999*x[i-1] + rnorm(1, sd = 0.01)
#   expect_lte(batchSize(x), 2)
# })

