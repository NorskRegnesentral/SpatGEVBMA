\name{gev.logscore}
\alias{gev.logscore}
\title{Compute the Log Score}


\description{This takes a large sample from a predictive distribution
  and computes the log score for a vector of observations, judging the
  fit of predictive distributions.
}
\usage{
gev.logscore(Y.obs, Y.samp)
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
 A scalar giving the mean log score over the observations
}
\author{
  Alex Lenkoski <alex@nr.no>
}
