\name{f.prime}
\alias{f.prime}

\title{First derivative of the posterior of a spatial GEV model with
  respect to a random effect in the location parameter.
}
\description{
  This internal function gives the first derivative of a spatial GEV
  model with respect to a random effect parameter in order to made a
  focused update of that variate.
}
\usage{
f.prime(tau, tau.hat, varsigma, xi, kappa, R)
}
\arguments{
  \item{tau}{
    Current value of the random effect
}
  \item{tau.hat}{
    Conditional mean of that random effect given the Gaussian process
  and the other random effects
}
  \item{varsigma}{
    Conditional variance of the random effect given the Gaussian process
}
  \item{xi}{
    Current shape parameter for that site
}
  \item{kappa}{
    Current precision parameter for that site
}
  \item{R}{
    Current vector of differences between the observations and the
  linear part of the location parameter at the site.
}
}
\value{
 A scalar giving the first derivative
}

\author{
  Alex Lenkoski <alex@nr.no>
}
