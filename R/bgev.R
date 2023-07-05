# PDF
pbgev = function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  #https://github.com/dcastrocamilo/bGEV/blob/master/Code/bGEVcode.R
  fix_lengths(x, mu, sigma, xi, p_a, p_b, s)
  g = get_gumbel_par(mu, sigma, xi, p_a, p_b)
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  p = pbeta((x - a) / (b - a), s, s)
  pgev(x, mu, sigma, xi) ^ p * pgev(x, g$mu, g$sigma, 0) ^ (1 - p)
}

# Quantile function
qbgev = function(p, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  fix_lengths(p, mu, sigma, xi, p_a, p_b, s)
  res = rep(NA, length(p))
  gumbel = which(p <= p_a)
  frechet = which(p >= p_b)
  mixing = which(p_a < p & p < p_b)
  if (any(gumbel)) {
    g = get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = qgev(p[gumbel], g$mu, g$sigma, 0)
  }
  if (any(frechet)) res[frechet] = qgev(p[frechet], mu[frechet], sigma[frechet], xi[frechet])
  if (any(mixing)) {
    res[mixing] = qbgev_mixing(p[mixing], mu[mixing], sigma[mixing], xi[mixing],
                               p_a[mixing], p_b[mixing], s[mixing])
  }
  res
}

# Random generation
rbgev = function(n, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  lengths = sapply(list(mu, sigma, xi), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qbgev(runif(n), mu, sigma, xi, p_a, p_b, s)
}

# Density using the usual parametrisation
dbgev = function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  fix_lengths(x, mu, sigma, xi, p_a, p_b, s)
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  res = rep(NA, length(x))
  gumbel = which(x <= a)
  frechet = which(x >= b)
  mixing = which(a < x & x < b)
  if (any(gumbel)) {
    g = get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = dgev(x[gumbel], g$mu, g$sigma, 0, log = log)
  }
  if (any(frechet)) res[frechet] = dgev(x[frechet], mu[frechet], sigma[frechet], xi[frechet], log = log)
  if (any(mixing)) {
    res[mixing] = dbgev_mixing(x[mixing], mu[mixing], sigma[mixing], xi[mixing],
                               p_a[mixing], p_b[mixing], s[mixing], log = log)
  }
  res
}

# Density using the new parametrisation
dbgev2 = function(x, q, sb, xi, alpha = 0.5, beta = 0.5, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  fix_lengths(x, q, sb, xi, p_a, p_b, s)
  tmp     = new.to.old(c(q,sb,xi), alpha = alpha, beta = beta)
  mu      = tmp$mu
  sigma   = tmp$sigma 
  a       = qgev(p_a, mu, sigma, xi)
  b       = qgev(p_b, mu, sigma, xi)
  res     = rep(NA, length(x))
  gumbel  = which(x <= a)
  frechet = which(x >= b)
  mixing  = which(a < x & x < b)
  if (any(gumbel)) {
    g           = get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = dgev(x[gumbel], g$mu, g$sigma, 0, log = log)
  }
  if (any(frechet)) res[frechet] = dgev(x[frechet], mu[frechet], sigma[frechet], xi[frechet], log = log)
  if (any(mixing)) {
    res[mixing] = dbgev_mixing(x[mixing], mu[mixing], sigma[mixing], xi[mixing],
                               p_a[mixing], p_b[mixing], s[mixing], log = log)
  }
  res
}

# Return levels
return_level_bgev = function(period, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  if (any(period <= 1)) warning("invalid period")
  p = ifelse(period > 1, 1 - 1 / period, NA)
  qbgev(p, mu, sigma, xi, p_a, p_b, s)
}

# Tool function for dgev
dbgev_mixing = function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  g = get_gumbel_par(mu, sigma, xi, p_a, p_b)
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  if (any(x <= a | x >= b)) stop("x is outside the domain for mixing")
  p = pbeta((x - a) / (b - a), s, s)
  p_der = dbeta((x - a) / (b - a), s, s) / (b - a)
  term1 = - p_der * (1 + xi * (x - mu) / sigma) ^ (-1 / xi)
  term2 = p / sigma * (1 + xi * (x - mu) / sigma) ^ (-1 / xi - 1)
  term3 = p_der * exp(- (x - g$mu) / g$sigma)
  term4 = (1 - p) / g$sigma * exp(- (x - g$mu) / g$sigma)
  term0 = p * log(pgev(x, mu, sigma, xi)) + (1 - p) * log(pgev(x, g$mu, g$sigma, 0))
  res = term0 + log(term1 + term2 + term3 + term4)
  if (!log) res = exp(res)
  res
}

# Tool function for qgev
qbgev_mixing = function(p, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, lower = 0, upper = 100) {
  if (any(p <= p_a | p >= p_b)) stop("p is outside the domain for mixing")
  res = vector("numeric", length(p))
  for (i in seq_along(p)) {
    f = function(x) (pbgev(x, mu, sigma, xi, p_a, p_b, s) - p)[i]
    sol = uniroot(f, lower = lower, upper = upper, extendInt = "upX")
    res[i] = sol$root
  }
  res
}



# Tool function to get parameters for G (Gumbel)
get_gumbel_par = function(mu, sigma, xi, p_a = .1, p_b = .2) {
  if (any(xi < 0)) stop("xi must be nonnegative")
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  sigma2 = (b - a) / log(log(p_a) / log(p_b))
  mu2 = a + sigma2 * log(-log(p_a))
  list(mu = mu2, sigma = sigma2)
}










#################################################
## Utility functions to work with the bGEV-GEV ##
#################################################
# By Silius M.V. and Daniela C.C.

new.to.old = function(par, alpha = 0.5, beta = 0.5){
  q = par[1] 
  s = par[2]
  xi = par[3]
  if(xi == 0){
    ell1 = log(-log(alpha))
    ell2 = log(-log(beta/2))
    ell3 = log(-log(1-beta/2))
  }else{
    ell1 = (-log(alpha))^(-xi)
    ell2 = (-log(beta/2))^(-xi)
    ell3 = (-log(1-beta/2))^(-xi)
  }
  if(xi == 0){
    sigma = s/(ell2 - ell3)
    mu    = q + sigma*ell1
  }else{
    mu    = q-s*(ell1 - 1)/(ell3-ell2)
    sigma = xi*s/(ell3-ell2)
  }
  list(mu = mu, sigma = sigma, xi = xi)
}
old.to.new = function(par, alpha = 0.5, beta = 0.5){
  mu    = par[1]
  sigma = par[2]
  xi    = par[3]
  
  qalpha = qgev(alpha, mu, sigma, xi)
  qbeta1 = qgev(beta/2, mu, sigma, xi)
  qbeta2 = qgev((1-beta/2), mu, sigma, xi)
  list(q = qalpha, s = qbeta2 - qbeta1, xi = xi)
}

fix_lengths = function(...) {
  call = match.call()
  varnames = sapply(call[-1], as.character)
  e = parent.frame()
  vars = lapply(varnames, get, envir = e)
  lengths = sapply(vars, length)
  max_length = max(lengths)
  if (any(max_length %% lengths != 0)) stop("Bad input lengths")
  for (i in seq_along(vars)) {
    if (lengths[i] < max_length) {
      assign(varnames[i], rep(vars[[i]], max_length / lengths[i]), envir = e)
    }
  }
  0
}



#############################################################################
## Negative log-likelihood functions associates to the GEV and bGEV models ##
#############################################################################
# By Daniela C.C.

# Neg log-lik of GEV using classical parametrisation
nllik.gev = function(par, x, log = TRUE){
  mu    = par[1]
  sigma = par[2]
  xi    = par[3]
  ll    = rep(NA, length(x))
  for(i in 1:length(x))
    ll[i] = dgev(x[i], mu, sigma, xi, log = log)
  -sum(ll)
}

# Neg log-lik of bGEV using classical parametrisation
nllik.bgev = function(par, x, p_a, p_b, s, log = TRUE){
  mu    = par[1]
  sigma = par[2]
  xi    = par[3]
  ll    = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dbgev(x[i], mu, sigma, xi, p_a = p_a, p_b = p_b, s = s, log = log)
    return(-sum(ll))
  }
}

# Neg log-lik of GEV using classical parametrisation with covariate-dependent location
nllik.gevx = function(par, x, w, log = TRUE){
  mu0   = par[1]
  mu1   = par[2]
  sigma = par[3]
  xi    = par[4]
  mu    = mu0 + mu1*w
  ll    = rep(NA, length(x))
  for(i in 1:length(x))
    ll[i] = dgev(x[i], mu[i], sigma, xi, log = log)
  
  -sum(ll)
}

# Neg log-lik of bGEV using classical parametrisation with covariate-dependent location
nllik.bgevx = function(par, x, w, p_a, p_b, s, log = TRUE){
  mu0   = par[1]
  mu1   = par[2]
  sigma = par[3]
  xi    = par[4]
  mu    = mu0 + mu1*w
  ll    = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dbgev(x[i], mu[i], sigma, xi, p_a = p_a, p_b = p_b, s = s, log = log)
    return(-sum(ll))
  }
}

# Neg log-lik of GEV using new parametrisation 
nllik.gev2 = function(par, x, alpha = 0.5, beta = 0.5, log = TRUE){
  q  = par[1]
  s  = par[2]
  xi = par[3]
  ll = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dgev2(x[i], q = q, s = s, xi = xi, alpha = alpha, beta = beta, log = log)
    return(-sum(ll))
  }
}

# Neg   log-lik of bGEV using new parametrisation
nllik.bgev2 = function(par, x, alpha = 0.5, beta = 0.5, p_a, p_b, s, log = TRUE){
  tmp   = new.to.old(par, alpha = alpha, beta = beta)
  mu    = tmp$mu
  sigma = tmp$sigma 
  ll    = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dbgev(x[i], mu, sigma, xi, p_a = p_a, p_b = p_b, s = s, log = log)
    return(-sum(ll))
  }
}


#################################################
## Functions to work with the GEV distribution ##
#################################################
# By Silius M.V. and Daniela C.C.

# PDF
pgev = function(x, mu, sigma, xi) {
  fix_lengths(x, mu, sigma, xi)
  ifelse(xi == 0,
         exp(-exp(- (x - mu) / sigma)),
         exp(-pmax(0, 1 + xi * (x - mu) / sigma) ^ (-1 / xi)))
}

# Quantile function
qgev = function(p, mu, sigma, xi) {
  fix_lengths(p, mu, sigma, xi)
  ifelse(xi == 0,
         mu - sigma * log(-log(p)),
         mu - sigma * (1 / xi) * (1 - (- log(p)) ^ (-xi)))
}

# Random generation
rgev = function(n, mu, sigma, xi) {
  lengths = sapply(list(mu, sigma, xi), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qgev(runif(n), mu, sigma, xi)
}

# Density using the usual parametrisation
dgev = function(x, mu, sigma, xi, log = FALSE) {
  fix_lengths(x, mu, sigma, xi)
  res = ifelse(xi == 0,
               -exp(- (x - mu) / sigma),
               -pmax(0, 1 + xi * (x - mu) / sigma) ^ (-1 / xi))
  res = res - log(sigma) +
    ifelse(xi == 0,
           - (x - mu) / sigma,
           ifelse(1 + xi * (x - mu) / sigma > 0,
                  - (1 / xi + 1) * log(1 + xi * (x - mu) / sigma),
                  -Inf))
  if (!log) res = exp(res)
  res
}

# Density using the new parametrisation
dgev2 = function(x, q, sb, xi, alpha = 0.5, beta = 0.5, log = FALSE) {
  tmp   = new.to.old(c(q,sb,xi), alpha = alpha, beta = beta)
  mu    = tmp$mu
  sigma = tmp$sigma 
  
  res = ifelse(xi == 0,
               -exp(- (x - mu) / sigma),
               -pmax(0, 1 + xi * (x - mu) / sigma) ^ (-1 / xi))
  res = res - log(sigma) +
    ifelse(xi == 0,
           - (x - mu) / sigma,
           ifelse(1 + xi * (x - mu) / sigma > 0,
                  - (1 / xi + 1) * log(1 + xi * (x - mu) / sigma),
                  -Inf))
  if (!log) res = exp(res)
  res
}

# Return levels
return_level_gev = function(period, mu, sigma, xi) {
  if (any(period <= 1)) warning("invalid period")
  p = ifelse(period > 1, 1 - 1 / period, NA)
  qgev(p, mu, sigma, xi)
}


#Made by TR
optim_gev=function(paramvec,y){
  loc=paramvec[1]
  scale=paramvec[2]
  shape=paramvec[3]
  
  liksum=0
  for(j in 1:length(y)){
    liksum=liksum+log(evd::dgev(y[j],loc=loc,scale=scale,shape=shape)+0.001)
  }
  
  return(-liksum)
}



gringorten <- function(x) {
  rank <- rank(x, na.last = "keep", ties.method = "first")
  len <- sum(!is.na(x))
  xx <- (rank - 0.44)/(len + 0.12)
  return(xx)
}


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
