\name{g.double.prime}
\alias{g.double.prime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The second derivative of a GEV distribution with respect to a
  random effect parameter on the precision kappa
}
\description{This function returns the second derivative of a GEV
  distribution with respect to the random effect tau_s^kappa.
  It is used in forming the proposal for the Metropolis-Hastings update
  of this parameter.
}
\usage{
g.double.prime(tau, tau.hat, varsigma, xi, kappa.hat, eps)
}

\arguments{
  \item{tau}{
    The current value of the random effect
  }
  \item{tau.hat}{
    The conditional fitted value of the random effect given
    the others according to the Gaussian Process
  }
  \item{varsigma}{
    The conditional variance of the fitted value according to the
    Gaussian Process
  }
  \item{xi}{
    The current value of the shape parameter
}
  \item{kappa.hat}{
    This is the linear bit of the precision for this location
}
  \item{eps}{
    The vector of residuals for this site, i.e. Y_ts - mu_s
}
}

\value{
  A scalar giving the second derivative
}
\author{
Alex Lenkoski <alex@nr.no>
}

