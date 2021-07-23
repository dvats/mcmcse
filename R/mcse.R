####################################################################
## Calculates the multivariate BM and MSV estimator for the covariance marix
## chain = MCMC output, matrix of size n x p
##
## method = which method to use. 
####################################################################
mcseR <- function(x, method = c("bm", "wbm", "bartlett", "tukey"), size = "sqroot", g = NULL, alpha = 0.90, warn = FALSE)
{ 
  chain <- as.matrix(x)
  
  if (is.function(g))
    chain <- t(apply(x, 1, g))

  ## Setting dimensions on the mcmc output. 
  n = dim(chain)[1]
  p = dim(chain)[2]
  method = match.arg(method)

  ## Initializing b_n 
  if(size == "sqroot")
  {
    b = floor(sqrt(n))
  } 
  else if(size == "cuberoot") {
    b = floor(n^(1/3))
  }
  else {
    if (!is.numeric(size) || size < 1 || size == Inf) 
        stop("'size' must be a finite numeric quantity larger than 1.")
    b = floor(size)
  }
  a = floor(n/b)

  if (a < p) {
      if (warn) 
          warning("Too few samples. Estimate need not be positive definite")
      if (n < 10) 
          return(NA)
  }
  
  
  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)
  
  ## Setting matrix sizes to avoid dynamic memory 
  sig.sum = matrix(0, nrow = p, ncol = p)
  
  ## only does bm and obm
  if(method != "bm" && method != "bartlett" && method != "tukey")
  {
    stop("No such method available")
  }
  
  ## Batch Means
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
    c <- exp(log(p) + log(a-1) - log(n) - log(a-p) + log(qf(alpha, p, a-p)))
  }


  if(method == "wbm")
  {
    a = floor(n/b)
    y.bar <- matrix(0,nrow = a, ncol = p)
    y.bar <- apply(chain, 2, function(x) sapply(1:a, function(k) mean(x[((k-1)*b+1):(k*b)])))
    for(i in 1:a)
    {
      sig.sum = sig.sum + tcrossprod(y.bar[i,] - mu.hat)
    }
 
    sig.mat1 <- b*sig.sum/(a-1)

    b <- b/2
    a = floor(n/b)
    y.bar <- matrix(0,nrow = a, ncol = p)
    y.bar <- apply(chain, 2, function(x) sapply(1:a, function(k) mean(x[((k-1)*b+1):(k*b)])))
    for(i in 1:a)
    {
      sig.sum = sig.sum + tcrossprod(y.bar[i,] - mu.hat)
    }
 
    sig.mat <- 2*sig.mat1 -  b*sig.sum/(a-1)

    c <- exp(log(p) + log(a-1) - log(n) - log(a-p) + log(qf(alpha, p, a-p)))
  }
  
  ## Modified Bartlett Window

  if(method == "bartlett")
  {
   chain <- scale(chain, center = mu.hat, scale = FALSE)
    dummy <- lapply(1:(b-1), function(j) (1 - j/b)*(  t(chain[1:(n-j), ])%*%chain[(j+1):n, ] 
                       + t(chain[(1+j):(n), ])%*%chain[1:(n-j), ]  ) )

   sig.sum <- t(chain)%*%chain + Reduce('+', dummy)

   sig.mat <- sig.sum/n
   c <- exp(log(p) + log(n-b) - log(n) - log(n - b - p +1) + log(qf(alpha, p, n - b - p + 1)))
  }

  if(method == "tukey")
  {
   chain <- scale(chain, center = mu.hat, scale = FALSE)
   dummy <- lapply(1:(b-1), function(j) ((.5)*(1 + cos(pi * j/b))*t(chain[1:(n-j), ]))%*%chain[(j+1):n, ] 
                      + ((.5)*(1 + cos(pi * j/b))*t(chain[(1+j):(n), ]))%*%chain[1:(n-j), ])

   sig.sum <- t(chain[1:(n), ])%*%chain[(1):n, ] + Reduce('+', dummy)
   sig.mat <- sig.sum/n
   c <- exp(log(p) + log(n-b) - log(n) - log(n - b - p +1) + log(qf(alpha, p, n - b - p + 1)))
  }


 # # c <- qchisq(alpha, p)/n
 #  c <- exp(log(p) + log(a-1) - log(n) - log(a-p) + log(qf(alpha, p, a-p)))
 
  dummy <- log(2) + (p/2)*log(pi*c) - log(p) - lgamma(p/2)
  log.det.sig <- sum(log(eigen(sig.mat, only.values = TRUE)$values))
  volume.sig <- exp(dummy + (1/2)*log.det.sig)
  volumes <- (volume.sig)^(1/p)
  return(list("cov" = sig.mat, "nsim" = n, "vol" = volumes))
}






