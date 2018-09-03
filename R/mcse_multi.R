# sourceCpp("mbmC.cpp")
# sourceCpp("msveC.cpp")

mcse.multi <- function(x, method = "bm", size = "sqroot", g = NULL, level = 0.95, large = FALSE)
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

  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)
  m <- 0
  ## Setting matrix sizes to avoid dynamic memory 
  sig.sum = matrix(0, nrow = p, ncol = p)
  
  if(method != "bm" && method != "wbm" && method != "lug" && method != "bartlett" && method != "tukey")
  {
    stop("No such method available")
  }
  
  ## Batch Means
  if(method == "bm")
  {
    sig.mat <- mbmC(chain, b)
    m <- a - 1
  }
  if(method == "wbm")
  {
    sig.mat <- 2*mbmC(chain, b) - mbmC(chain, floor(b/2))
    m <- a - 1
  }
  if(method == "lug")
  {
    sig.mat <- 2*mbmC(chain, b) - mbmC(chain, floor(b/3))
    m <- a - 1
  }
  ## Modified Bartlett Window

  if(method == "bartlett")
  {
   chain <- scale(chain, center = mu.hat, scale = FALSE)
   sig.mat <- msveC(chain, b, "bartlett")
   m <- n - b
  }

  ## tukey window
  if(method == "tukey")
  {
   chain <- scale(chain, center = mu.hat, scale = FALSE)
   sig.mat <- msveC(chain, b, "tukey")
   m <- n - b
  }


  if(m - p +1 <=0)
  {
    warning("Not enough samples. Estimate is not positive definite. ")
    pth.vol <- NaN
  } else
  {
    sig.eigen <- eigen(sig.mat, only.values = TRUE)$values
    if(min(sig.eigen) <= 0)
    {
      warning("You either need more samples or x is not full column rank")
    }
    log.dethalf.pth<- (1/(2*p))*sum(log(sig.eigen))
    crit <- ifelse(large, qchisq(level, df = p)/n,
              exp(log(p) + log(m) - log(n) - log(m-p+1) + log(qf(level, p, m-p+1))) )
    foo <- (1/p)*log(2) + (1/2)*log(pi*crit) - (1/p)*log(p) - (1/p)*lgamma(p/2)
    pth.vol <- exp(foo+log.dethalf.pth)
}
  return(list("cov" = sig.mat, "vol" = pth.vol, "est" = mu.hat, "nsim" = n, 
      "method" = method, "large" = large, "size" = b))
}






