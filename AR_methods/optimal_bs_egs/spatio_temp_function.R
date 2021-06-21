library(spBayes)
data("NETemp.dat")
ne.temp <- NETemp.dat


spatio_temp <- function(N = 5e4, verbose = FALSE)
{
##take a chunk of New England
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]##subset first 2 years (Jan 2000 - Dec. 2002)
y.t <- ne.temp[,4:15]
N.t <- ncol(y.t) ##number of months
n <- nrow(y.t) ##number of observation per months

coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))

##set starting and priors
p <- 2 #number of regression parameters in each month
starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
"sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
"sigma.eta"=diag(rep(0.01, p)))

tuning <- list("phi"=rep(6, N.t))

priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
"phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
"sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
"tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
"sigma.eta.IW"=list(2, diag(0.001,p)))

##make symbolic model formula statement for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)
n.samples <- N

m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
starting=starting, tuning=tuning, priors=priors, get.fitted =FALSE,
cov.model="exponential", n.samples=n.samples, n.report=1e4, verbose = verbose)


chain <- cbind(m.1$p.beta.0.samples,
m.1$p.beta.samples,
m.1$p.theta.samples,
m.1$p.sigma.eta.samples[ ,-3],
t(m.1$p.u.samples))

foo <- rep(1:N.t, each = dim(y.t)[1])
foo <- paste("t",foo, "s", rep(1:dim(y.t)[1], N.t), sep = "")
foo <- paste("u.", foo, sep = "")
colnames(chain)[(p + N.t*p + 3*N.t + p*(p+1)/2 + 1):dim(chain)[2]] <- foo

return(chain)
}


