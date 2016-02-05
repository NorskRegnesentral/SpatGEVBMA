\name{l.double.prime}
\alias{l.double.prime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The second derivative of a Gaussian process with respect to the
  parameter lambda.

}
\description{This utility returns the value of the second derivative of
  the posterior of lambda in a Gaussian process.  It is used to make a
  focused M-H proposal when updating the lambda parameter.
}

\usage{
l.double.prime(tau, alpha, lambda, D, a, b)
}

\arguments{
  \item{tau}{
    The vector of current random effects.
}
  \item{alpha}{
    The current state of the precision of the Gaussian process
}
  \item{lambda}{
    The current state of lambda
}
  \item{D}{
    The distance matrix of the relevant entries in tau
}
  \item{a}{
    The prior parameter for lambda 
}
  \item{b}{
    The second prior parameter for lambda
}
}

\value{

  This returns a scalar that is then used in the focused update.  Its a
  pretty internal function.
}
\author{
  Alex Lenkoski <alex@nr.no>
}
