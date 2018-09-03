\name{estvssamp}
\alias{estvssamp}
\title{Create a plot that shows how Monte Carlo estimates change with increasing sample size.}
\usage{
  estvssamp(x, g = mean, main = "Estimates vs Sample Size",
    add = FALSE, ...)
}
\arguments{
  \item{x}{a sample vector.}

  \item{g}{a function such that \eqn{E(g(x))} is the
  quantity of interest. The default is \code{g =
  \link{mean}}.}

  \item{main}{an overall title for the plot. The default is
  \dQuote{\code{Estimates vs Sample Size}}.}

  \item{add}{logical. If \code{TRUE}, add to a current
  plot.}

  \item{\dots}{additional arguments to the plotting
  function.}
}
\value{
  \code{NULL}
}
\description{
  Create a plot that shows how Monte Carlo estimates change
  with increasing sample size.
}
\examples{
\dontrun{
estvssamp(x, main = expression(E(beta)))
estvssamp(y, add = TRUE, lty = 2, col = "red")}
}

