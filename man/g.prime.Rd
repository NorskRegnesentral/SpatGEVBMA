\name{g.prime}
\alias{g.prime}

\title{The first derivative of the posterior density of a spatial GEV
  model with respect to a given random effect on the precision parameter.

}
\description{This returns the first derivative of the posterior density
  of a spatial GEV model with respect to a random effect on the
  precision parameter.   It is used in forming the Metropolis-Hastings
  update of this parameter.
}
\usage{
g.prime(tau, tau.hat, varsigma, xi, kappa.hat, eps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{
    Current value of the random effect
}
  \item{tau.hat}{
    The conditional mean of the random effect given the others and the
    current Gaussian process parameters.
  }
  \item{varsigma}{
    The conditional variance of the random effect based on the Gaussian
    process parameters
}
  \item{xi}{
    The current shape parameter for this location
}
  \item{kappa.hat}{
    The linear part of the precision parameters
}
  \item{eps}{
    The vector of residuals based on the observations at this site and
    the associated location parmeter
}
}

\value{
  A scalar giving the first derivative, which is used to form the
  Metropolis-Hasting update.
}

\author{
 Alex Lenkoski <alex@nr.no>
}
