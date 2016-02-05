\name{spatial.gev.bma}
\alias{spatial.gev.bma}

\title{Run an MCMC to fit a hierarchical spatial generalized extreme
  value (GEV) model with the option for Bayesian model averaging (BMA)

}
\description{This is the main function of the spatial.gev.bma package.
  It runs an MCMC to sample the posterior distribution of a spatial GEV
  model that takes a number of covariates for linear interactions as
  well as a spatial component to model dependent overdispersion.
}
\usage{
spatial.gev.bma(Y.list, X.all, S, n.reps, prior.user = NULL,
full = FALSE, fixed.xi = NULL, print.every = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.list}{
    This is a list of length n, where n is the total number of sites you
    have data for.  Each element of the list is a vector of observations
    at that site.
}
  \item{X.all}{
    This is an n by p matrix where n is the total number of sites and p
    is the total number of variables under consideration.  Please supply
    your own constant in the first column.  Note that this model
    implicitly assumes that at a given site, the covariates are
    constant throughout the observation period.  Email me if you need
    something more flexible.
}
  \item{S}{
    This is an n by 2 matrix with the coordinates of each of your n
    sites.
}
  \item{n.reps}{
   The total number of repititions (i.e. includes burn in) that you
   would like to run the MCMC.
}
  \item{prior.user}{
    A list containing any priors that the user would like to specify.
}
  \item{full}{
    A boolean indicating if you'd like to force all variables to be
    included (TRUE) or whether you want to perform BMA (FALSE)
}
  \item{fixed.xi}{
    If set to NULL, the shape parameter is estimated like the location
    and scale.  If a numeric value is provided, the shape parameter is
    fixed at this level.
}
  \item{print.every}{
    A simple tool to tell you how many iterations you've run.  If you
    set it to 0 then no progress reports will be given. If it is a
    positive number, say x, then every x iterations it will tell you the
    total number of iterations.
}
}

\value{
  \item{TAU}{An n by reps by 3 array of the states of the random effects
  at each iteration}
  \item{THETA}{An reps by p by 3 array of the states of the linear terms
    at each iteration}

  \item{ALPHA}{A reps by 3 matrix of states of the alpha term at each
    iteration}
  \item{LAMBDA}{A reps by 3 matrix of states of the lambda term at each
    iteration}
}


\author{
  Alex Lenkoski <alex@nr.no>
}

\examples{
data(norway)
attach(norway)
##To replicate our results, change 2 to 2e5 below
a <- spatial.gev.bma(Y.list,X,S,2)  
}
