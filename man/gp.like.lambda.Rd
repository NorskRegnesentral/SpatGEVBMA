\name{gp.like.lambda}
\alias{gp.like.lambda}

\title{The likelihood of a Gaussian process used to initialize the
  lambda parameter

}
\description{This is a utility that is used to initialize the lambda
  parameter in a Gaussian process when initializing the entire spatial
  GEV model.  Its rather single service.
}
\usage{
gp.like.lambda(lambda, alpha, tau, D)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
    The lambda parameter
}
  \item{alpha}{
    The alpha parmater of a gaussian process
}
  \item{tau}{
    The vector of random effects
}
  \item{D}{
    The distance matrix for the entries in tau
}
}
\value{
 The log density of the Gaussian process.

 }