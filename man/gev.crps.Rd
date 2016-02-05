\name{gev.crps}
\alias{gev.crps}
\title{Compute the Continuous Rank Probability Score (CRPS)

}
\description{This takes a large sample from a predictive distribution
  and computes the CRPS for a number of observations.  It is used in
  judging the fit of predictive distributions.

}
\usage{
gev.crps(Y.obs, Y.samp)
}

\arguments{
  \item{Y.obs}{
    A vector of observed values
}
  \item{Y.samp}{
    A large sample from the posterior predictive distribution.
}
}

\value{
 A scalar giving the CRPS
}
\author{
  Alex Lenkoski <alex@nr.no>
}
