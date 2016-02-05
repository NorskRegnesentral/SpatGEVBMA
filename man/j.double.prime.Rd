\name{j.double.prime}
\alias{j.double.prime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The second derivative of a spatial GEV with respect to a random
  effect in the shape parameter

}
\description{This internal function helps with the MCMC update of a RE
  on the shape parameter

}
\usage{
j.double.prime(tau, tau.hat, varsigma, kappa, xi.hat, eps)
}
\arguments{
  \item{tau}{
    Current value of the RE
}
  \item{tau.hat}{
    Conditional mean given the GP
}
  \item{varsigma}{
    Conditional variance given the GP
}
  \item{kappa}{
    Current precision for this site
}
  \item{xi.hat}{
    Linear part of the shape parameter
}
  \item{eps}{
    The residuals between the observation and the location parameter at
    this site
}
}

\value{
  A scalar giving the second derivative
}

\author{
  Alex Lenkoski <alex@nr.no>
}
