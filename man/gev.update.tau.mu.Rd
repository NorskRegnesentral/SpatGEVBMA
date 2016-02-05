\name{gev.update.tau.mu}
\alias{gev.update.tau.mu}

\title{Internal function to update the random effects of the location
  parameter in a Spatial GEV model.

}
\description{This internal function takes the state of a Spatial GEV
  model, and updates the random effects for the location parameter
}
\usage{
gev.update.tau.mu(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{An object holding the current state information for a Spatial
    GEV model}
}

\value{
  An updated state for use in the Spatial GEV model.
}
\author{
  Alex Lenkoski <alex@nr.no>
}
