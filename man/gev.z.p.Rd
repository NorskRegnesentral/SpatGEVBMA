\name{gev.z.p}
\alias{gev.z.p}
\title{Calculate the 1/p return level for a GEV distribution
}
\description{This simple function returns the 1/p return level for a GEV distribution

}
\usage{
gev.z.p(p, mu, sigma, xi)
}

\arguments{
  \item{p}{The level on a 0 to 1 scale.  Note that this is associated
  with the 1/p return level.
}
  \item{mu}{
    The location parameter of a GEV distribution
}
  \item{sigma}{
    The scale parameter of a GEV distribution
}
  \item{xi}{
    The shape parameter of a GEV distribution
}
}

\value{
 A scalar giving the associated return level.
}
\author{
Alex Lenkoski <alex@nr.no>
}
