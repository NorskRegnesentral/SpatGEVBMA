\name{gev.like}
\alias{gev.like}
\title{The log likelihood of a GEV distribution

}
\description{Returns the log likelihood of a GEV distribution accoding
  to the parameterization used in spatial.gev.bma

}
\usage{
gev.like(Y, mu, kappa, xi)
}
\arguments{
  \item{Y}{The outcome at which the log likelihood should be evaluated
}
  \item{mu}{The location parameter
}
  \item{kappa}{The precision parameter
}
  \item{xi}{The shape parameter
}
}
\value{
A scalar giving the log-likelihood
}
\author{
  Alex Lenkoski <alex@nr.no>
}
