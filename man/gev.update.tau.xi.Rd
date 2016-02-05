\name{gev.update.tau.xi}
\alias{gev.update.tau.xi}
\title{Update the random effects for the shape parameter in a spatial
  GEV model
}
\description{
  This goes through the random effects in the shape parameter of a
  hierarchical spatial GEV model and updates them individually using
  focused Metropolis-Hastings proposals.
}
\usage{
gev.update.tau.xi(G)
}

\arguments{
  \item{G}{
    The current state of the model
}
}

\value{
  A new state of the model where the random effects have been updated.
}

\author{
  Alex Lenkoski <alex@nr.no>
}
