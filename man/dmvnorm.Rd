\name{dmvnorm}
\alias{dmvnorm}
\title{
Log density of a multivariate normal distribution
}
\description{
This returns the important bits of the log of a multivariate normal
distribution, for use in Metropolis-Hastings ratios.  One of these
functions you wish they would just include in base so you don't have to
keep rolling your own or calling some other package.
}
\usage{
dmvnorm(x, mu, Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A length p vector at which you would like the density to be evaluated
}
  \item{mu}{A length p vector which is the mean of the MV normal
}
  \item{Sigma}{A p x p matrix for the covariance of the MV normal
}
}
\value{
The log density of a multivariate normal
}
\author{
Alex Lenkoski <alex@nr.no>
}
