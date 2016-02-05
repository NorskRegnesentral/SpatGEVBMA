\name{gev.process.results}
\alias{gev.process.results}

\title{Outputs some tables from the results of Spatial GEV MCMC run

}
\description{This take the output from a spatial GEV run and returns
  some handy tables given posterior probabilities, effective sample
  sizes, acceptance probabilities, etc.

}
\usage{
gev.process.results(R, burn = 100)
}
\arguments{
  \item{R}{
    The output from the spatial.gev.bma function
}
  \item{burn}{
    An optional burn-in amount to be discarded before forming results.
}
}

\value{
  Returns a number of tables.
}

\author{
  Alex <alex@nr.no>
}
