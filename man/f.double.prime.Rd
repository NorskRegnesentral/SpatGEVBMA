\name{f.double.prime}
\alias{f.double.prime}
\title{Second derivative of the posterior distribution of a spatial GEV
  with respect a location random effect
}
\description{
  This is the second derivative of the conditional posterior probability
  of a sptial GEV model with respect to a given random effect in the
  location parameter.  It is used to form the proposal for a
  Metropolis-Hastings update of the given parameter.
}
\usage{
f.double.prime(tau, tau.hat, varsigma, xi, kappa, R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{
    The current value of the given random effect
}
  \item{tau.hat}{
    The conditional mean of the random effect given the other random
    effects and the Gaussian process parameters
}
  \item{varsigma}{
    The conditional variance of the random effect given the Gaussian
    process parameters
}
  \item{xi}{
    The current shape parameter at this location
}
  \item{kappa}{
    The current precision parameter at this location
}
  \item{R}{
    The residual between the observed values and the linear part of the
    location parameter
}
}
\value{
  A scalar that gives the value of the derivative.  This is then used to
  form the relevant proposal.
}
\author{
  Alex Lenkoski <alex@nr.no>
}
