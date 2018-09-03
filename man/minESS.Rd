\name{minESS}
\alias{minESS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum effective sample size required for stable estimation as described in Vats et al. (2015).
}
\description{
The function calculates the minimum effective sample size required for a specified relative tolerance level. This function can also calculate the relative precision in estimation for a given estimated effective sample size.
}
\usage{
minESS(p, alpha = .05, eps = .05, ess = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{p}{dimension of the estimation problem.}
\item{alpha}{confidence level}
\item{eps}{tolerance level. The \code{eps} value is ignored is \code{ess} is not \code{NULL}}
\item{ess}{Estimated effective sample size. Usually the output value from \code{multiESS}.}
}
\details{
The minimum effective samples required when estimating a vector of length p, with 100(1-\eqn{\alpha})\% confidence and tolerance of \eqn{\epsilon} is

\deqn{\mbox{mESS} \geq \frac{2^{2/p} \pi}{(p \Gamma(p/2))^{2/p}} \frac{\chi^2_{1-\alpha, p}}{\epsilon^2} }

The above equality can also be used to get \eqn{\epsilon} from an already obtained estimate of mESS.
}
\value{
  By default function returns the minimum effective sample required for a given \code{eps} tolerance. If \code{ess} is specified, then the value returned is the \code{eps} corresponding to that \code{ess}.
}


\examples{
minESS(p = 5)
}

\references{
Gong, L., and Flegal, J. M. A practical sequential stopping rule for high-dimensional Markov chain Monte Carlo. Journal of Computational and Graphical Statistics (to appear).

Vats, D., Flegal, J. M., and, Jones, G. L Multivariate Output Analysis for Markov chain Monte Carlo, \emph{arXiv preprint arXiv:1512.07713 (2015)}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\seealso{
  \code{\link{multiESS}}, which calculates multivariate effective sample size using a Markov chain and a function g.

  \code{\link{ess}} which calculates univariate effective sample size using a Markov chain and a function g.
}