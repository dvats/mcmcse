
#' Compute Monte Carlo standard errors for expectations. 
#'
#' @param x a vector of values from a Markov chain of length n.
#' @param size represents the batch size in \dQuote{\code{bm}} and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code
#' {tukey}}. Default is \code{NULL} which implies that an optimal batch size is calculated using the 
#'   \code{batchSize} function. Can take character values of \dQuote{\code{sqroot}} and \dQuote{\code{cuberoot}} or any numeric
#'   value between 1 and n/2. \dQuote{\code{sqroot}} means size is floor(n^(1/2)) and \dQuote{\code{cuberoot}} means size is
#'   floor(n^(1/3)).
#' @param g a function such that \eqn{E(g(x))} is the quantity of interest. The default is
#'   \code{NULL}, which causes the identity function to be used.
#' @param method any of \dQuote{\code{bm}},\dQuote{\code{obm}},\dQuote{\code{bartlett}}, \dQuote{\code{tukey}}. \dQuote{\code{bm}}
#' represents batch means estimator, \dQuote{\code{obm}} represents overlapping batch means estimator with, \dQuote{\code{bartlett}}
#' and \dQuote{\code{tukey}} represents the modified-Bartlett window and the Tukey-Hanning windows for spectral variance estimators.
#'   
#' @param warn a logical value indicating whether the function should issue a warning if the sample
#'   size is too small (less than 1,000).
#' @param r The lugsail parameters (\code{r}) that converts a lag window into its lugsail
#'   equivalent. Larger values of \code{r} will typically imply less underestimation of \dQuote{\code{cov}},
#'   but higher variability of the estimator. Default is \code{r = 3} and \code{r = 1,2} are
#'   good choices. \code{r > 5} is not recommended.
#' 
#' @return \code{mcse} returns a list with three elements:
#'         \item{est}{an estimate of \eqn{E(g(x))}.}
#'         \item{se}{the Monte Carlo standard error.}
#'         \item{nsim}{The number of samples in the input Markov chain.}
#'         
#' @references
#' Flegal, J. M. (2012) Applicability of subsampling bootstrap methods in Markov chain Monte Carlo.
#' In Wozniakowski, H. and Plaskota, L., editors, \emph{Monte Carlo and Quasi-Monte Carlo Methods
#' 2010} (to appear). Springer-Verlag.
#'
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov
#' chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070.
#'
#' Flegal, J. M. and Jones, G. L. (2011) Implementing Markov chain Monte Carlo: Estimating with
#' confidence. In Brooks, S., Gelman, A., Jones, G. L., and Meng, X., editors, \emph{Handbook of
#' Markov Chain Monte Carlo}, pages 175--197. Chapman & Hall/CRC Press.
#'
#' Flegal, J. M., Jones, G. L., and Neath, R. (2012) Markov chain Monte Carlo estimation of
#' quantiles. \emph{University of California, Riverside, Technical Report}.
#'
#' Jones, G. L., Haran, M., Caffo, B. S. and Neath, R. (2006) Fixed-width output analysis for Markov
#' chain Monte Carlo. \emph{Journal of the American Statistical Association}, \bold{101}, 1537--154.
#' 
#' @seealso
#' \code{\link{mcse.mat}}, which applies \code{mcse} to each column of a matrix or data frame.
#' 
#' \code{\link{mcse.multi}}, for a multivariate estimate of the Monte Carlo standard error.
#'
#' \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard errors for quantiles.
#' 
#' @examples
#'
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' out = BVN_Gibbs(n, mu, sigma)
#' x = out[,1]
#' mcse(x)
#' mcse.q(x, 0.1)
#' mcse.q(x, 0.9)
#'
#' # Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using overlapping batch means.
#'
#' mcse(x, method = "obm")
#' mcse.q(x, 0.1, method = "obm")
#' mcse.q(x, 0.9, method = "obm")
#'
#' # Estimate E(x^2) with MCSE using spectral methods.
#'
#' g = function(x) { x^2 }
#' mcse(x, g = g, method = "tukey")
#'
#' @export

mcse <- function(x, size = NULL, g = NULL, r=3, method = c("bm", "obm", "bartlett", "tukey"), warn = FALSE)
{
    x <- as.numeric(x)
    n = length(x)
    method = match.arg(method)
    calls <- mcse.multi(x, method = method, r=r, size = size, g = g, adjust = FALSE, blather = FALSE)
    se <- as.numeric(sqrt(calls$cov / n))
    mu.hat <- calls$est
    value = list("est" = mu.hat, "se" = se, "nsim" = n)
    value
}

#' Apply \code{mcse} to each column of a matrix or data frame of MCMC samples.
#'
#' @param x a matrix of values from a Markov chain of size n x p.
#' @param size represents the batch size in \dQuote{\code{bm}} and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code
#' {tukey}}. Default is \code{NULL} which implies that an optimal batch size is calculated using the 
#'   \code{batchSize} function. Can take character values of \dQuote{\code{sqroot}} and \dQuote{\code{cuberoot}} or any numeric
#'   value between 1 and n/2. \dQuote{\code{sqroot}} means size is floor(n^(1/2)) and \dQuote{\code{cuberoot}} means size is
#'   floor(n^(1/3)).
#' @param g a function such that \eqn{E(g(x))} is the quantity of interest. The default is
#'   \code{NULL}, which causes the identity function to be used.
#' @param method any of \dQuote{\code{bm}},\dQuote{\code{obm}},\dQuote{\code{bartlett}}, \dQuote{\code{tukey}}. \dQuote{\code{bm}}
#' represents batch means estimator, \dQuote{\code{obm}} represents overlapping batch means estimator with, \dQuote{\code{bartlett}}
#' and \dQuote{\code{tukey}} represents the modified-Bartlett window and the Tukey-Hanning windows for spectral variance estimators.
#'   
#' @param r The lugsail parameters (\code{r}) that converts a lag window into its lugsail
#'   equivalent. Larger values of \code{r} will typically imply less underestimation of \dQuote{\code{cov}},
#'   but higher variability of the estimator. Default is \code{r = 3} and \code{r = 1,2} are
#'   good choices. \code{r > 5} is not recommended.
#' 
#' @return \code{mcse.mat} returns a matrix with \code{ncol(x)} rows and two columns. The row names
#'   of the matrix are the same as the column names of \code{x}. The column names of the matrix are
#'   \dQuote{\code{est}} and \dQuote{\code{se}}. The \eqn{j}th row of the matrix contains the result
#'   of applying \code{mcse} to the \eqn{j}th column of \code{x}.
#'   
#' @seealso
#' \code{\link{mcse}}, which acts on a vector.
#' 
#' \code{\link{mcse.multi}}, for a multivariate estimate of the Monte Carlo standard error.
#'
#' \code{\link{mcse.q}} and \code{\link{mcse.q.mat}}, which compute standard errors for quantiles.
#' 
#' @export

mcse.mat = function(x, size = NULL, g = NULL, method = c("bm", "obm", "bartlett", "tukey"), r=3)
{
    if (! is.matrix(x) && ! is.data.frame(x))
        stop("'x' must be a matrix or data frame.")
    num = ncol(x)
    vals = matrix(NA, num, 2)
    colnames(vals) = c("est", "se")
    rownames(vals) = colnames(x)
    res = apply(x, 2, mcse, size = size, g = g, method = method, r=r)
    for (i in 1:num)
        vals[i, ] = c(res[[i]]$est, res[[i]]$se)
    vals
}

quant = function(input, q) { quantile(input, prob = q, type = 1, names = FALSE) }

#' Compute Monte Carlo standard errors for quantiles.
#'
#' @param x a vector of values from a Markov chain.
#' @param q the quantile of interest.
#' @param size the batch size. The default value is \dQuote{\code{sqroot}}, which uses the square
#'   root of the sample size. A numeric value may be provided if \dQuote{\code{sqroot}} is not
#'   satisfactory.
#' @param g a function such that the \eqn{q}th quantile of the univariate distribution function of
#'   \eqn{g(x)} is the quantity of interest. The default is \code{NULL}, which causes the identity
#'   function to be used.
#' @param method the method used to compute the standard error. This is one of \dQuote{\code{bm}}
#'   (batch means, the default), \dQuote{\code{obm}} (overlapping batch means), or
#'   \dQuote{\code{sub}} (subsampling bootstrap).
#' @param warn a logical value indicating whether the function should issue a warning if the sample
#'   size is too small (less than 1,000).
#'   
#' @return \code{mcse.q} returns a list with three elements:
#'         \item{est}{an estimate of the \eqn{q}th quantile of the univariate distribution function of \eqn{g(x)}.}
#'         \item{se}{the Monte Carlo standard error.}
#'         \item{nsim}{The number of samples in the input Markov chain.}
#'         
#' @references
#' Flegal, J. M. (2012) Applicability of subsampling bootstrap methods in Markov chain Monte Carlo.
#' In Wozniakowski, H. and Plaskota, L., editors, \emph{Monte Carlo and Quasi-Monte Carlo Methods
#' 2010} (to appear). Springer-Verlag.
#'
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov
#' chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070.
#'
#' Flegal, J. M. and Jones, G. L. (2011) Implementing Markov chain Monte Carlo: Estimating with
#' confidence. In Brooks, S., Gelman, A., Jones, G. L., and Meng, X., editors, \emph{Handbook of
#' Markov Chain Monte Carlo}, pages 175--197. Chapman & Hall/CRC Press.
#'
#' Flegal, J. M., Jones, G. L., and Neath, R. (2012) Markov chain Monte Carlo estimation of
#' quantiles. \emph{University of California, Riverside, Technical Report}.
#'
#' Jones, G. L., Haran, M., Caffo, B. S. and Neath, R. (2006) Fixed-width output analysis for Markov
#' chain Monte Carlo. \emph{Journal of the American Statistical Association}, \bold{101}, 1537--154
#' .
#' 
#' @seealso
#' \code{\link{mcse.q.mat}}, which applies \code{mcse.q} to each column of a matrix or data frame.
#'
#' \code{\link{mcse}} and \code{\link{mcse.mat}}, which compute standard errors for expectations.
#' @examples
#'
#' ## Bivariate Normal with mean (mu1, mu2) and covariance sigma
#' n <- 1e3
#' mu = c(2, 50)
#' sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' out = BVN_Gibbs(n, mu, sigma)
#' x = out[,1]
#'
#' # Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using batch means.
#'
#' mcse(x)
#' mcse.q(x, 0.1)
#' mcse.q(x, 0.9)
#'
#' # Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using overlapping batch means.
#'
#' mcse(x, method = "obm")
#' mcse.q(x, 0.1, method = "obm")
#' mcse.q(x, 0.9, method = "obm")
#'
#' # Estimate E(x^2) with MCSE using spectral methods.
#'
#' g = function(x) { x^2 }
#' mcse(x, g = g, method = "tukey")
#'
#' @export

mcse.q <- function(x, q, size = NULL, g = NULL, method = c("bm", "obm", "sub"), warn = FALSE) {
    
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
    
    if (! is.numeric(q) || q <= 0 || q >= 1)    {
        stop("'q' must be from (0, 1)")
    }
    
    xi.hat = quant(x, q)
    
    if(method == "bm")  {
        var.hat = mcseqbm(x, b, xi.hat)
        f.hat.junk = density(x, from = xi.hat, to = xi.hat, n = 1)
        f.hat = f.hat.junk$y
        se = sqrt(var.hat / n) / f.hat
    } else if(method == "obm")  {
        var.hat = mcseqobm(x, b, xi.hat)
        f.hat.junk = density(x, from = xi.hat, to = xi.hat, n = 1)
        f.hat = f.hat.junk$y
        se = sqrt(var.hat / n) / f.hat
    } else  {
        var.hat = mcseqsub(x, b, q, quant)
        se = sqrt(var.hat / n)
        
    }
    
    
    value = list("est" = xi.hat, "se" = se, "nsim" = n)
    value
}

#' Apply \code{mcse.q} to each column of a matrix or data frame of MCMC samples.
#'
#' @param x a matrix or data frame with each row being a draw from the multivariate distribution of
#'   interest.
#' @param q the quantile of interest.
#' @param size the batch size. The default value is \dQuote{\code{sqroot}}, which uses the square
#'   root of the sample size. \dQuote{\code{cuberoot}} will cause the function to use the cube root
#'   of the sample size. A numeric value may be provided if \dQuote{\code{sqroot}} is not
#'   satisfactory.
#' @param g a function such that the \eqn{q}th quantile of the univariate distribution function of
#'   \eqn{g(x)} is the quantity of interest. The default is \code{NULL}, which causes the identity
#'   function to be used.
#' @param method the method used to compute the standard error. This is one of \dQuote{\code{bm}}
#'   (batch means, the default), \dQuote{\code{obm}} (overlapping batch means), or
#'   \dQuote{\code{sub}} (subsampling bootstrap).
#'   
#' @return \code{mcse.q.mat} returns a matrix with \code{ncol(x)} rows and two columns. The row
#'   names of the matrix are the same as the column names of \code{x}. The column names of the
#'   matrix are \dQuote{\code{est}} and \dQuote{\code{se}}. The \eqn{j}th row of the matrix contains
#'   the result of applying \code{mcse.q} to the \eqn{j}th column of \code{x}.
#'   
#' @seealso \code{\link{mcse.q}}, which acts on a vector.
#'
#' \code{\link{mcse}} and \code{\link{mcse.mat}}, which compute standard errors for expectations.
#' @export

mcse.q.mat = function(x, q, size = NULL, g = NULL, method = c("bm", "obm", "sub"))
{
    if (! is.matrix(x) && ! is.data.frame(x))
        stop("'x' must be a matrix or data frame.")
    num = ncol(x)
    vals = matrix(NA, num, 2)
    colnames(vals) = c("est", "se")
    rownames(vals) = colnames(x)
    res = apply(x, 2, mcse.q, q = q, size = size, g = g, method = method)
    for (i in 1:num)
        vals[i, ] = c(res[[i]]$est, res[[i]]$se)
    vals
}

#' Create a plot that shows how Monte Carlo estimates change with increasing sample size.
#'
#' @param x a sample vector.
#' @param g a function such that \eqn{E(g(x))} is the quantity of interest. The default is \code{g =
#'   \link{mean}}.
#' @param main an overall title for the plot. The default is \dQuote{\code{Estimates vs Sample
#'   Size}}.
#' @param add logical. If \code{TRUE}, add to a current plot.
#' @param \dots additional arguments to the plotting function.
#' 
#' @return \code{NULL}
#' 
#' @examples
#' \dontrun{
#' estvssamp(x, main = expression(E(beta)))
#' estvssamp(y, add = TRUE, lty = 2, col = "red")}
#' @export

estvssamp = function(x, g = mean, main = "Estimates vs Sample Size", add = FALSE,...)
{
    if (length(x) < 100)
        size = 1
    else
        size = length(x) %/% 100
    n = seq(size, length(x), by = size)
    est = c()
    for (j in n)
        est = c(est, g(x[1:j]))
    if (add)
        lines(n, est,...)
    else
        plot(n, est, main = main, type = "l", xlab = "Sample Size", ylab = "MC Estimate",...)
}

#' Univariate estimate effective sample size (ESS) as described in Gong and Flgal (2015).
#'
#'Estimate effective sample size (ESS) as described in Gong and Flegal (2015).
#'
#' @details 
#' ESS is the size of an iid sample with the same variance as the current sample for estimating the expectation of g. ESS is given by
#' \deqn{ESS = n \frac{\lambda^{2}}{\sigma^{2}}} where \eqn{\lambda^{2}} is the sample variance and
#' \eqn{\sigma^{2}} is an estimate of the variance in the Markov chain central limit theorem. This is by default
#' a batch means estimator, but the default can be changed with the `method` argument.
#'
#' @param x a matrix or data frame of Markov chain output. Number of rows is the Monte
#'   Carlo sample size.
#' @param ... arguments passed on to the mcse.mat function. For example method = \dQuote{\code{tukey}} and size =
#'   \dQuote{\code{cuberoot}} can be used.
#' @param g a function that represents features of interest. \code{g} is applied to each row of x and thus
#'  \code{g} should take a vector input only. Ifcode{g} is \code{NULL}, \code{g} is set to be identity, which is estimation
#'   of the mean of the target density.
#'   
#' @return The function returns the estimated effective sample size for each component of \code{g}.
#' 
#' @references
#' Gong, L. and Flegal, J. M. (2015) A practical sequential stopping rule for high-dimensional
#' Markov chain Monte Carlo, Journal of Computational and Graphical Statistics.
#' 
#' @seealso 
#' \code{\link{minESS}}, which calculates the minimum effective samples required for the problem.
#' \code{\link{multiESS}}, which calculates multivariate effective sample size using a Markov chain
#' and a function \code{g}.
#' 
#' @export

ess <- function(x, g = NULL, ...)
{
    chain <- as.matrix(x)
    
    if (is.function(g)) 
        chain <- t(apply(x, 1, g))
    
    n <- dim(chain)[1]
    p <- dim(chain)[2]
    
    lambda <- apply(chain, 2, var)
    
    sigma <- as.numeric((mcse.mat(chain,...)[,2])^2*n)
    
    return(n*lambda/sigma)
}

