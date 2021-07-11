library(Rcpp)
sourceCpp('../batchsize.cpp')

###################################################
### Final update of batchSize using thresholding
### on aR coefficients and option to use only the
### tail of the chain for acf calculation
###################################################
batchSize_final <- function(x, method = "bm", g = NULL, last_size = 5e4, fast = TRUE) {
  
  if(!is.numeric(x))
    stop("'x' must be numeric") # the chain must be numeric
  if(any(is.na(x)))
    stop("NAs found")     # stop in NAs found
  
  p <- ncol(x)
  n <- sum(!is.na(x[,1])) # number of non-missing rows
  chain <- as.matrix(x)
  
  if(!is.matrix(chain) && !is.data.fradim(me(chain)))
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
  
  order.max <- min(p, n - 1L, floor(10 * log10(n))) # Maximum order up to which AR is fit
  xacf = matrix(, nrow = order.max+1, ncol = p)

  if(fast)  {                                       # Use only the tail of the chain to calculate acf
    last = min(n, last_size)
    chain2 = chain[(n-last+1):n,]
    xacf = sapply(1:p, function(i) acf(chain2[,i], type = "covariance", lag.max = order.max, plot = FALSE,
                                       demean=TRUE, na.action = na.pass)$acf)
  }

  else  {                                           # use the entire chain for acf calculation
    xacf = sapply(1:p, function(i) acf(chain[,i], type = "covariance", lag.max = order.max, plot = FALSE,
                                       demean=TRUE, na.action = na.pass)$acf)
  }

  threshold = qnorm((1.95)/2)/sqrt(n)              # threshold used in confidence interaval calculation
  b = batchsize_cpp(n, p, xacf, order.max, method, threshold)
  return(b)
}
