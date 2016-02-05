\name{gev.update.tau.kappa}
\alias{gev.update.tau.kappa}

\title{Update the random effects of the precision parameter in a spatial
  GEV model

}
\description{This goes through and updates, via M-H, each component in
  the random effects of a spatial GEV model

}
\usage{
gev.update.tau.kappa(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
    The current state of a spatial GEV model
}
}

\value{
  A new state of the spatial GEV model where the random effects have
  been updated
}
\author{
  Alex Lenkoski <alex@nr.no>
}
