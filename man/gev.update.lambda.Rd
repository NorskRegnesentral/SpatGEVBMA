\name{gev.update.lambda}
\alias{gev.update.lambda}

\title{Update the lambda parameter in a Gaussian Process

}
\description{
  This constructs a focused proposal to make a metropolis-hasting update
  of the lambda parameter in a Gaussian process.  It is agnostic to
  which model component in the spatial GEV the GP belongs to.
}
\usage{
gev.update.lambda(tau, alpha, lambda, D, a, b)
}
\arguments{
  \item{tau}{
    A vector of the random effects
}
  \item{alpha}{
    The global precision of the Gaussian process
}
  \item{lambda}{
    The current value of lambda
}
  \item{D}{
    The matrix of distances between the points associated with the
    vector tau
}
  \item{a}{
    The parameter for the prior of lambda
}
  \item{b}{
    The second parameter for the prior of lambda
}
}

\value{
 A scalar which is the updated value of lambda.
}

\author{
  Alex Lenkoski <alex@nr.no>
}
