\name{gev.update.theta}
\alias{gev.update.theta}

\title{Update the linear parameters in a spatial GEV model

}
\description{This internal function updates the linear regression terms
  of a given component of a spatial GEV model.  Due to blocking, the
  method is agnostic to which component it is updating

}
\usage{
gev.update.theta(Y, X, M, alpha, lambda, D, beta.0, Omega.0)
}

\arguments{
  \item{Y}{The dependent valuesn

}
  \item{X}{The matrix of covariates

}
  \item{M}{The current model for this component

}
  \item{alpha}{The overall variance of the GP

}
  \item{lambda}{The scale of the GP

}
  \item{D}{The distance matrix of the GP

}
  \item{beta.0}{The prior mean for the covariates

}
  \item{Omega.0}{The prior precision for the covariates

}
}

\value{
  Returns a new set of linear regression terms
}

\author{
  Alex Lenkoski <alex@nr.no>
}
