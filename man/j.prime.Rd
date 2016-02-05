\name{j.prime}
\alias{j.prime}

\title{The first derivative of the posterior density of a spatial GEV
  model with respect to a random effect parameter on the shape.
}
\description{This internal function calculates the first derivative of a
  spatial GEV model with respect to a random effect parameter on the
  shape and is used during the Metropolis-Hastings update of this
  parameter conditional on everything else.

}
\usage{
j.prime(tau, tau.hat, varsigma, kappa, xi.hat, eps)
}
\arguments{
  \item{tau}{
  Current value of the random effect
}
  \item{tau.hat}{
    Conditional mean of the random effect given all others and the
    Gaussian process parameters.
}
  \item{varsigma}{
    Conditional variance of the random effect according to the Gaussian process
}
  \item{kappa}{
    Current value for the GEV scale at this site.
}
  \item{xi.hat}{
    Current linear part of the GEV shape at this site.
}
  \item{eps}{
    Vector of residuals for observations of this site given the current
    location value.
}
}

\value{
  A scalar giving the first derivative of this density
}

\author{
 Alex Lenkoski <alex@nr.no>
}
