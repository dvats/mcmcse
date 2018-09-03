\name{ess}
\alias{ess}
\title{Estimate effective sample size (ESS) as described in Gong and Felgal (2015).}
\usage{
  ess(x, g = NULL, ...)
}
\arguments{
\item{x}{a matrix or data frame of Markov chain output. Number of rows is the Monte Carlo sample size.}
  \item{...}{arguments passed on to the \code{mcse.mat} function. For example \code{method = "tukey"} and \code{size = "cuberoot"} can be used.}
  % \item{imse}{logical. If \code{TRUE}, use an approach that
  % is analogous to Geyer's initial monotone positive
  % sequence estimator (IMSE), where correlations beyond a
  % certain lag are removed to reduce noise.}

  % \item{verbose}{logical. If \code{TRUE} and \code{imse =
  % TRUE}, inform about the lag at which truncation occurs,
  % and warn if the lag is probably too small.}
\item{g}{a function that represents features of interest. g is applied to each row of \code{x} and thus g should take a vector input only. If g is \code{NULL}, g is set to be identity, which is estimation of the mean of the target density.}
}
\value{
  The function returns the estimated effective sample size.
}
\description{
  Estimate effective sample size (ESS) as described in Gong and Flegal (2015).
}
\details{
  ESS is the size of an iid sample with the same variance
  as the current sample. ESS is given by
  \deqn{\mbox{ESS}=n \frac{\lambda^2}{\sigma^2},}{ESS = n  \lambda^2/\sigma^2,} where
\eqn{\lambda^2} is the sample variance and \eqn{\sigma^2} is an estimate of the variance in the CLT. This is by default the batch means estimator, but the default can be changed with the \code{method} argument.
}
\references{

Gong, L. and Flegal, J. M. (2015) A practical sequential stopping rule for high-dimensional Markov chain Monte Carlo \emph{Journal of Computational and Graphical Statistics}.

  % Kass, R. E., Carlin, B. P., Gelman, A., and Neal, R.
  % (1998) Markov chain Monte Carlo in practice: A roundtable
  % discussion. \emph{The American Statistician}, \bold{52},
  % 93--100.

  % Robert, C. P. and Casella, G. (2004) \emph{Monte Carlo
  % Statistical Methods}. New York: Springer.

  % Geyer, C. J. (1992) Practical Markov chain Monte Carlo.
  % \emph{Statistical Science}, \bold{7}, 473--483.
}

\seealso{
  \code{\link{minESS}}, which calculates the minimum effective samples required for the problem.
  
  \code{\link{multiESS}}, which calculates multivariate effective sample size using a Markov chain and a function g.
}