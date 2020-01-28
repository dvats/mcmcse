# sourceCpp("inseq.cpp")

mcse.initseq <- function(x, g = NULL, adjust = FALSE, blather = FALSE)
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
  
  ## Setting dimensions on the mcmc output
  n = dim(chain)[1]
  p = dim(chain)[2]
  
  ## Overall means of the mcmc output
  mu.hat <- colMeans(chain)

  ## Initial Sequence Estimator(s)
  res <- inseq(chain, adjust)
  
  ##sig=initial sequence estimator without asjustment
  sig <- res$Sig
  
  ##sig.adj=initial sequence estimator with asjustment, if adjust=T
  ##       =NULL, if adjust=F
  if(adjust)
  {
    sig.adj <- res$Sigadj
  }else
  {
    sig.adj <- NULL
  }

  
  if(blather)
  {
    return(list("cov" = sig, "cov.adj"=sig.adj,
              "est" = mu.hat, "nsim" = n, "adjust" = adjust)) 
  } else{
    return(list("cov" = sig, "cov.adj"=sig.adj,
              "est" = mu.hat)) 
  }

}



