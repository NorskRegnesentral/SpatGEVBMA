\name{l.prime}
\alias{l.prime}

\title{First derivative of a GP with respect to lambda

}
\description{This internal function gives the first derivative of a GP
  with respect to the parameter lambda
}

\usage{
  l.prime(tau, alpha, lambda, D, a, b)
}
\arguments{
  \item{tau}{
    The current vector of REs
}
  \item{alpha}{
    The current overall variance of the process
}
  \item{lambda}{
    The current value of lambda
}
  \item{D}{
    The distance matrix
}
  \item{a}{
    The first prior parameter for lambda
}
  \item{b}{
    The second prior parameter for lambda
}
}

\value{
A scalar giving the derivative
}

\author{
  Alex Lenkoski <alex@nr.no>
}
