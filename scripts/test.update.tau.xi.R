rm(list = ls())

# Run the following line after writing new cpp-code to test.
Rcpp::compileAttributes("~/jullum/Prosjekter/SpatGEV/Git_SpatialGEVBMA/SpatialGEVBMA")

#Rcpp::compileAttributes("~/jullum/Prosjekter/SpatGEV/old")


# Run the following in the terminal:
# R CMD INSTALL ~/Prosjekter/SpatGEV/Git_SpatialGEVBMA/SpatialGEVBMA

####setwd("H:\Prosjekter\SpatGEV\Git_SpatialGEVBMA")  # For windows


library("SpatialGEVBMA")

j.double.prime <- function(tau, tau.hat, varsigma, kappa, xi.hat, eps)# The old function here
{
  ## What an ugly function
  h <- 1 + kappa * eps * (xi.hat + tau)
  if(any(h < 0))return(-Inf)
  f.1 <- (xi.hat + tau + 1)/(xi.hat + tau) * log(h)
  f.2 <-  exp(-(xi.hat + tau)^(-1) * log(h))
  
  f.1.dot <- -log(h)/(xi.hat + tau)^2 + (xi.hat + tau + 1)/(xi.hat + tau) * h^(-1) * eps * kappa
  f.2.dot <- f.2 * ( log(h)/(xi.hat + tau)^2 - h^(-1) * eps * kappa / (xi.hat + tau) )
  g.1.dot <- -2 * (xi.hat + tau)^(-3) * log(h) + (xi.hat + tau)^(-2) * h^(-1) * kappa * eps
  g.2.dot <- -h^(-1) * eps * kappa * (xi.hat + tau)^(-2) - (xi.hat + tau + 1)/(xi.hat + tau) * h^(-2) * eps^2 * kappa^2
  g.3.dot.1 <- f.2.dot * (log(h) * (xi.hat + tau)^(-2) )
  g.3.dot.2 <- f.2 * ( -2 * log(h) * (xi.hat + tau)^(-3) + h^(-1) * kappa * eps * (xi.hat + tau)^(-2) )
  g.3.dot <- g.3.dot.1 + g.3.dot.2
  g.4.dot.1 <- f.2.dot * (h^(-1) * kappa * eps * (xi.hat + tau)^(-1))
  g.4.dot.2 <- -f.2 * eps * kappa * ( h^(-1) * (xi.hat + tau)^(-2) + h^(-2) * eps * kappa * (xi.hat + tau)^(-1) )
  g.4.dot <- g.4.dot.1 + g.4.dot.2
  res <- sum(g.1.dot) - sum(g.2.dot) - sum(g.3.dot) + sum(g.4.dot) - 1/varsigma
  return(res)
}

j.prime <- function(tau, tau.hat, varsigma, kappa, xi.hat, eps)
{
  h <- 1 + kappa * eps * (xi.hat + tau)
  if(any(h < 0))return(-Inf)
  f.1 <- (xi.hat + tau + 1)/(xi.hat + tau) * log(h)
  f.1.dot <- -log(h)/(xi.hat + tau)^2 + (xi.hat + tau + 1)/(xi.hat + tau) * h^(-1) * eps * kappa
  f.2 <-  exp(-(xi.hat + tau)^(-1) * log(h))
  f.2.dot <- f.2 * (log(h)/(xi.hat + tau)^2 - h^(-1) * eps * kappa / (xi.hat + tau) )
  res <- -sum(f.1.dot) - sum(f.2.dot) - (tau - tau.hat)/varsigma
  return(res)
}


data(norway)
attach(norway)

set.seed(1)
G <- gev.init(Y.list, X, S, NULL, FALSE, NULL)
for(i in 1:5)G <- gev.update(G)



##------ Unload --------
n.s <- G$n.s
theta.xi <- G$theta.xi

C <- 1/G$alpha.xi * exp(-1/G$lambda.xi * G$D)
diag(C) <- diag(C) + 1e-5
C.inv <- solve(C)
tau <- G$tau.xi
accept <- 0
G$accept.tau.xi <- rep(0,n.s)

s<- 1 # Test value
  Y.s <- G$Y.list[[s]]
  mu.s <- sum(G$theta.mu * G$X[s,]) + G$tau.mu[s]
  eps.s <- Y.s - mu.s
  kappa.s <- sum(G$theta.kappa * G$X[s,]) + G$tau.kappa[s]
  xi.hat <- sum(G$theta.xi * G$X[s,])
  
  tau.hat.s <- sum(-C.inv[s,-s]/C.inv[s,s] * tau[-s])
  varsigma.s <- 1/C.inv[s,s] ##precision matrix stuff
  

#j.s <- j.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)



##--------------------
# Some kind of testing right here.



jj.s <- j.double.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
jj.s.arma <- j_double_prime_new(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)

j.s <- j.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
j.s.arma <- j_prime_new(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)

jj.s
jj.s.arma
### Identical

j.s
j.s.arma
## Identical

kappa.s=-100
h <- 1 + kappa.s * eps.s * (xi.hat + tau[2])

jj.s <- j.double.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
jj.s.arma <- j_double_prime_new(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)

j.s <- j.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
j.s.arma <- j_prime_new(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)

jj.s
jj.s.arma
### Identical

j.s
j.s.arma


### Timing:

kappa.s <- sum(G$theta.kappa * G$X[s,]) + G$tau.kappa[s]# Reset kappa value


ptm=proc.time()
s=1
for (i in 1:10^5)
{
  jj.s <- j.double.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
}
proc.time()-ptm
#user  system elapsed 
#7.532   0.000   7.544

ptm=proc.time()
s=1
for (i in 1:10^5)
{
  jj.s <- j_double_prime_new(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
}
proc.time()-ptm
# user  system elapsed 
# 2.364   0.000   2.377 


ptm=proc.time()
s=1
for (i in 1:10^5)
{
  j.s <- j.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
}
proc.time()-ptm
#user  system elapsed 
#2.700   0.000   2.715 

ptm=proc.time()
s=1
for (i in 1:10^5)
{
  j.s <- j_prime_new(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
}
proc.time()-ptm
#user  system elapsed 
#1.556   0.000   1.580


