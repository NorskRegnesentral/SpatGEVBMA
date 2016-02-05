\name{logdet}
\alias{logdet}
\title{Returns the log determinant for a symmetric positive definite matrix.

}
\description{
  This internal function returns the log determinant of a positive
  definite matrix by using the trace of a the Cholesky trick.  Its way
  faster than what you find in base.
 
}
\usage{
logdet(A)
}

\arguments{
  \item{A}{A symmetric positive definite matrix.

}
}

\value{
A scalar giving the log determinant of A
}
\author{
Alex Lenkoski <alex@nr.no>
}
