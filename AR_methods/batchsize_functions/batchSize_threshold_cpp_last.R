library(Rcpp)
sourceCpp('../batchsize.cpp')

batchSize_threshold_cpp_last <- function(x, method = "bm", g = NULL, last_size = 5e4) {
  
  if(!is.numeric(x))
    stop("'x' must be numeric")
  if(any(is.na(x)))
    stop("NAs found")
  p <- ncol(x)

  n <- sum(!is.na(x[,1])) # number of non-missing rows
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
  order.max <- min(p, n - 1L, floor(10 * log10(n))) 
  xacf = matrix(, nrow = order.max+1, ncol = p)
  last = min(n, last_size)
  chain2 = chain[(n-last+1):n,]
  xacf = sapply(1:p, function(i) acf(chain2[,i], type = "covariance", lag.max = order.max, plot = FALSE,
                                     demean=TRUE, na.action = na.pass)$acf)

  threshold = qnorm((1.95)/2)/sqrt(n)
  b = batchsize_cpp(n, p, xacf, order.max, method, threshold)
  return(b)
}