\name{gev.update.M}
\alias{gev.update.M}

\title{Sample a new model from the current model for any linear
  regression system

}
\description{
  This uses a conditional Bayes factor (CBF) to update a model in a
  linear system given a current model and other information in a spatial
  GEV model.  Note that it is agnostic to which part of the framework
  (location, precision, scale) you are updating.
}
\usage{
gev.update.M(Y, X, M, alpha, lambda, D, beta.0, Omega.0)
}

\arguments{
  \item{Y}{
    The current dependent variable, calculated relative to the linear
  plus random effect terms of the given component.
}
  \item{X}{
    The matrix of covariates
}
  \item{M}{
    The current model.  A subset of (1, ..., p) where p is the number of
  columns in X
}
  \item{alpha}{
  The precision term of the Gaussian process for this component of the model
}
  \item{lambda}{
    The length term of the Gaussian process for this component of the model
}
  \item{D}{
    The distance matrix used in the Gaussian process
}
  \item{beta.0}{
   The prior mean on the linear regression terms
}
  \item{Omega.0}{
    The prior covariance on the linear regression terms
}
}

\value{
  This returns an updated model, which is a vector that is a subset of
  (1, ..., p).
}
\author{
  Alex Lenkoski <alex@nr.no>
}
