library(Rcpp)
sourceCpp('../batchsize.cpp')

batchSize_threshold_cpp <- function(x, method = "bm", g = NULL) {
  
  if(!is.numeric(x))
    stop("'x' must be numeric")
  if(any(is.na(x) != is.na(x[,1]))) stop("NAs in 'x' must be the same row-wise")
  p <- ncol(x)
  # xm <- colMeans(x, na.rm=TRUE)
  # x <- sweep(x, 2L, xm, check.margin=FALSE)
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
  for(i in 1:p) {
    xacf[,i] = acf(chain[,i], type = "covariance", lag.max = order.max, plot = FALSE,
                   demean=TRUE, na.action = na.pass)$acf
  }
  # xacf <- acf(chain, type = "covariance", lag.max = order.max, plot = FALSE,
  #               demean=TRUE, na.action = na.pass)$acf
  ci = 0.95
  threshold = qnorm((1 + ci)/2)/sqrt(N)
  b = batchsize_cpp(n, p, xacf, order.max, method, threshold)
  return(b)
}