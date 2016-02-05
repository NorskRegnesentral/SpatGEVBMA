\name{gev.init}
\alias{gev.init}

\title{Initilizes a state object for a Spatial GEV distribution

}
\description{This utility function just initializes the state object
  that is then passed around during the course of an MCMC for a
  hierarchical spatial GEV distribution. It starts the object in a nice
  spot and fills in the relevant items.

}
\usage{
gev.init(Y.list, X.all, S, prior.user, full, fixed.xi)
}

\arguments{
  \item{Y.list}{The list of Y values

}
  \item{X.all}{The matrix of covariates

}
  \item{S}{The matrix of spatial locations

}
  \item{prior.user}{The user-specified prior parameters, if there are any

}
  \item{full}{A boolean indicating if BMA is to be performed

}
  \item{fixed.xi}{An indication if the parameter xi is to be estimated (NULL)
  or set to a given value.
}
}

\value{
  This returns a state object that is then used for the MCMC by the main
  updating routines.
}

\author{
 Alex Lenkoski <alex@nr.no>
}
