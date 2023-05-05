library(INLA)
library(evgam)
rm(list = ls())
library(data.table)
loc_thea = "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/"
setwd(loc_thea)
source("R/gev.R")
source("R/temporal_spatgev.R")

#Upload hourly data (change 60 to 1440 to get daily data):
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")

#Select some covariates:
amax_data=amax_data[,.(year,masl,stid,lon,lat,
                       wetterdays_JJA,wetterdays_annual,
                       precip_JJA,precip_annual,
                       temp_JJA,temp_annual,y)]

#Select locations in the Bergen area:
#amax_data=amax_data[lon<7&lat<62 & lat>58]
amax_data[,intercept:=1]

#---------------------------------------#
n = 1000
x = rnorm(n, sd=0.5) # we generate values for x from a N(0,0.5ˆ2) dist.
eta.x = 1 + 0.4*x
spread = 0.3
tail = 0.1
p.alpha = 0.5
p.beta = 0.25
par = giveme.gev.par(q = eta.x, sbeta = spread, alpha = p.alpha, beta = p.beta,
                     xi = tail)
y = numeric(n)
for(i in 1:n)
  y[i] = rgev(1, loc = par$mu[i], scale = par$sigma, shape = par$xi)

hyper.spread = list(initial = 1,
                    fixed=FALSE,
                    prior = "loggamma",
                    param = c(3, 3))
tail.interval = c(0, 0.5)
tail.intern = map.tail(tail, tail.interval, inverse=TRUE)

hyper.tail = list(initial = tail.intern,
                  prior = "pc.gevtail",
                  param = c(7, tail.interval),
                  fixed= FALSE)

hyper.tail = list(initial = if (tail == 0.0) -Inf else tail.intern,
                  prior = "pc.gevtail",
                  param = c(7, tail.interval),
                  fixed= if (tail == 0.0) TRUE else FALSE)

hyper.bgev = list(spread = hyper.spread,
                  tail = hyper.tail)

control.bgev = list(q.location = p.alpha,
                    q.spread = p.beta,
                    # quantile levels for the mixing part
                    q.mix= c(0.05, 0.20),
                    # the Beta(s1, s2) mixing distribution parameters.
                    # Hard-coded for the moment: s1=s2=5
                    beta.ab = 5)

null.matrix = matrix(nrow = n, ncol= 0)
spread.x = null.matrix
tail.x = null.matrix

data.bgev = data.frame(y = y, intercept = 1, x = x, spread.x = spread.x, tail.x = tail.x)
formula = inla.mdata(y, spread.x, tail.x) ~ -1 + intercept + x
#formula = inla.mdata(y, 1, 1) ~ -1 + intercept + x
formula = inla.mdata(y) ~ -1 + intercept + x

r1 = inla(formula,
          family = "bgev",
          data = data.bgev,
          control.family = list(hyper = hyper.bgev,
                                control.bgev = control.bgev),
          control.predictor = list(compute = TRUE),
          control.fixed = list(prec=1000),
          control.compute = list(cpo = TRUE),
          control.inla = list(int.strategy = "eb"),
          verbose=FALSE, safe=TRUE)

r1$summary.fixed
r1$summary.hyperpar
r1$summary.random


amax2= data.table(scale(amax_data[,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl)]))
amax2[,intercept:=1]

amax3=round(amax_data[,.(wetterdays_JJA,precip_JJA,lon,lat,temp_JJA,masl,y)],1)
amax3[,precip_JJA:=precip_JJA/1000]
amax3[,intercept:=1]

formula = inla.mdata(y,amax3[,.(lat,lat)] ,amax3[,.(lat,lat)]) ~ -1 + intercept + wetterdays_JJA+precip_JJA+temp_JJA+masl+lon+lat


r1 = inla(formula,
          family = "bgev",
          data = amax2,
          control.family = list(hyper = hyper.bgev,
                                control.bgev = control.bgev),
          control.predictor = list(compute = TRUE),
          control.fixed = list(prec=1000),
          control.compute = list(cpo = TRUE),
          control.inla = list(int.strategy = "eb"),
          verbose=FALSE, safe=TRUE)

linpred=r1$summary.fixed$`0.5quant`
hyperpar=r1$summary.hyperpar$`0.5quant`

r1$summary.fixed
r1$summary.hyperpar

par = giveme.gev.par(q = r1$summary.fitted.values$mean, sbeta = r1$summary.hyperpar$mean[1], alpha = p.alpha, beta = p.beta,
                     xi = tail)


#eta.x=

#par = giveme.gev.par(q = eta.x, sbeta = spread, alpha = p.alpha, beta = p.beta,
 #                    xi = tail)


#https://cran.r-project.org/web/packages/SpatialGEV/vignettes/SpatialGEV-vignette.html






library(evd)
giveme.gev.par = function(q, sbeta, alpha, beta, xi)
{
  .mu = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    c = -log(alpha)
    if (all(xi > 0.0)) {
      tmp0 = (c^(-xi) - 1)/xi
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(q - (sbeta/dbeta) * tmp0)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      tmp0 = log(c)
      return(q + (sbeta/dbeta) * tmp0)
    } else {
      stop("mixed case not implemented")
    }
  }
  .sigma = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    if (all(xi > 0.0)) {
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(sbeta/dbeta)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      return(sbeta/dbeta)
    } else {
      stop("mixed case not implemented")
    }
  }
  return(list(mu = .mu(q, sbeta, alpha, beta, xi),
              sigma = .sigma(q, sbeta, alpha, beta, xi),
              xi = xi))
}


map.tail = function(x, interval, inverse = FALSE) {
  if (!inverse) {
    return (interval[1] + (interval[2] - interval[1]) * exp(x)/(1.0 + exp(x)))
  } else {
    return (log((x-interval[1])/(interval[2]-x)))
  }
}

