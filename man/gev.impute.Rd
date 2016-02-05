\name{gev.impute}
\alias{gev.impute}
\title{Given the output of the MCMC, return a number of samples for a new site.
}
\description{
  This takes the object returned from spatial.gev.bma and new
  information for a given site and returns a large sample of values
  from the posterior predictive distribution for that site.
}
\usage{
gev.impute(R, X.drop, S.drop, burn = NULL, n.each = NULL)
}
\arguments{
  \item{R}{
    This is the object returned from spatial.gev.bma
}
  \item{X.drop}{
    This is the vector of covariates for the site that has been left out.
}
  \item{S.drop}{
    This is the 2 dimensional vector of locations for this site
}
  \item{burn}{
    Optional number of repetitions to drop from R as burn in.  Defaults
    to 10 percent.
}
  \item{n.each}{
    Optional number of samples from each repetition under
    consideration.  Defaults to a number such that around 1 million
    samples are returned.
}
}

\value{
  Returns a large number of sample from the posterior predictive
  distribution.
}
\author{
  Alex Lenkoski <alex@nr.no>
}
