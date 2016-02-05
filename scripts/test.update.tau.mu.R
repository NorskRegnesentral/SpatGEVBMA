
rm(list = ls())

# Run the following line after writing new cpp-code to test.
#Rcpp::compileAttributes("/nr/samba/user/jullum/Prosjekter/FlomQ/Git_SpatialGEVBMA/SpatialGEVBMA")

# Run the following in the terminal:
# cd /nr/samba/user/jullum/Prosjekter/FlomQ/Git_SpatialGEVBMA
# R CMD INSTALL --library="/nr/samba/user/jullum/Rlibs/new_packages" SpatialGEVBMA



library("SpatialGEVBMA", lib.loc="/nr/samba/user/jullum/Rlibs/new_packages")



data(norway)
attach(norway)

set.seed(1)
G <- gev.init(Y.list, X, S, NULL, FALSE, NULL,nonspatial=FALSE,log.kappa=FALSE)
for(i in 1:5)G <- gev.update(G)



##------ Unload --------
n.s <- G$n.s
theta.mu <- G$theta.mu

C <- exp(-1/G$lambda.mu * G$D) / G$alpha.mu
diag(C) <- diag(C) + 1e-5
C.inv <- solve(C)
tau <- G$tau.mu
G$accept.tau.mu <- rep(0, n.s)
s<-15

  Y.s <- G$Y.list[[s]]
  fit.s <- sum(theta.mu * G$X[s,])
  R.s <- Y.s - fit.s
  kappa.s <- sum(G$theta.kappa * G$X[s,]) + G$tau.kappa[s]
  xi.s <- sum(G$theta.xi * G$X[s,]) + G$tau.xi[s]
  tau.hat.s <- sum(-C.inv[s,-s]/C.inv[s,s] * tau[-s])
  varsigma.s <- 1/C.inv[s,s] ##precision matrix stuff
  


### Testing:

f.s <- f.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
ff.s <- f.double.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)


f.s.arma <- f_prime_new(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
ff.s.arma <- f_double_prime_new(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)


f.s
f.s.arma

ff.s
ff.s.arma


f.s <- f.prime(tau[s], tau.hat.s, varsigma.s, 0, kappa.s, R.s)
ff.s <- f.double.prime(tau[s], tau.hat.s, varsigma.s, 0, kappa.s, R.s)


f.s.arma <- f_prime_new(tau[s], tau.hat.s, varsigma.s, 0, kappa.s, R.s)
ff.s.arma <- f_double_prime_new(tau[s], tau.hat.s, varsigma.s, 0, kappa.s, R.s)


f.s
f.s.arma

ff.s
ff.s.arma


### Identical


ptm=proc.time()
for (i in 1:(5*10^5))
{
  f.s <- f.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
}
proc.time()-ptm
#user  system elapsed 
#4.832   0.008   4.849 


ptm=proc.time()
for (i in 1:(5*10^5))
{
  f.s.arma <- f_prime_new(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
}
proc.time()-ptm
#user  system elapsed 
#4.409   0.000   4.414 

# 10% speedup

ptm=proc.time()
for (i in 1:(5*10^5))
{
  ff.s <- f.double.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
}
proc.time()-ptm
#user  system elapsed 
#5.492   0.004   5.507 


ptm=proc.time()
for (i in 1:(5*10^5))
{
  ff.s.arma <- f_double_prime_new(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
}
proc.time()-ptm
#user  system elapsed 
#4.916   0.000   4.924 

# 10% speedup
