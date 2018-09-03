minESS <- function(p, alpha = .05, eps = .05, ess = NULL)
{
	crit <- qchisq(1-alpha, p)
	foo <- 2/p

	if(is.null(ess) == TRUE)
	{
		logminESS <- foo*log(2) + log(pi) - foo*log(p) - foo*lgamma(p/2) - 2*log(eps) + log(crit)
		names(logminESS) <- "minESS"
		return(round(exp(logminESS)))
	}else{
		if(is.numeric(ess) == FALSE)
		{
			stop("Only numeric entry allowed for ess")
		}
		logEPS <- .5*foo*log(2) + .5*log(pi) - .5*foo*log(p) - .5*foo*lgamma(p/2) - .5*log(ess) + .5*log(crit)
		names(logEPS) <- "Epsilon"
		return(exp(logEPS))
	}
}