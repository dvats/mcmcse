# sourceCpp("inseq.cpp")

mcse.initseq <- function(x, adjust = FALSE, g = NULL, level = 0.95)
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
  
  #calculate volume of the confidence region to the pth root
  crit <- qchisq(level,df=p)/n
  #log scale
  foo <- -(log(p/2) + lgamma(p/2))/p + log(pi*crit)/2
  det2p <- log(det(sig))/2/p
  ##vol=volume to the pth root without adjustment
  vol <- exp(foo + det2p)
  
  ##vol.adj=volume to the pth root with adjustment, if adjust=T
  ##       =NULL, if adjust=F
  if(adjust)
  {
    det2p.adj <- log(det(sig.adj))/2/p
    vol.adj <- exp(foo + det2p.adj)
  }else
  {
    vol.adj <- NULL
  }
  
  return(list("cov" = sig, "cov.adj"=sig.adj,
              "vol"=vol, "vol.adj"=vol.adj,
              "est" = mu.hat, "nsim" = n, "adjust" = adjust))
}