\name{make.D}
\alias{make.D}
\title{Form the distance matrix for use in a Gaussian Process

}
\description{
This takes two matrices containing the (lat,lon) locations of two
collections of points and returns the matrix of euclidean distances
between them.
}
\usage{
make.D(x.1, x.2)
}

\arguments{
  \item{x.1}{
    An n1 by 2 matrix of the locations of the first set of points
}
  \item{x.2}{
    An n2 by 2 matrix of the locations of the second set of points
}
}
\value{
  Returns an n1 by n2 matrix with the relevant Euclidean distances
}

\author{
  Alex Lenkoski <alex@nr.no>
}
