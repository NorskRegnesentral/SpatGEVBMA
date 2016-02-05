
rm(list = ls())

# Run the following line after writing new cpp-code to test.
#Rcpp::compileAttributes("/nr/samba/user/jullum/Prosjekter/FlomQ/Git_SpatialGEVBMA/SpatialGEVBMA")

# Run the following in the terminal:
# cd /nr/samba/user/jullum/Prosjekter/FlomQ/Git_SpatialGEVBMA
# R CMD INSTALL --library="/nr/samba/user/jullum/Rlibs/new_packages" SpatialGEVBMA



library("SpatialGEVBMA", lib.loc="/nr/samba/user/jullum/Rlibs/new_packages")


## Unload
data(norway)
attach(norway)

set.seed(1)
G <- gev.init(Y.list, X, S, NULL, FALSE, NULL,nonspatial=FALSE,log.kappa=FALSE)
for(i in 1:5)G <- gev.update(G)


s=3
n.s <- G$n.s
theta.kappa <- G$theta.kappa

C <- 1/G$alpha.kappa * exp(-1/G$lambda.kappa * G$D)
diag(C) <- diag(C) + 1e-5
C.inv <- solve(C)
tau <- G$tau.kappa



Y.s <- G$Y.list[[s]]
mu.s <- sum(G$theta.mu * G$X[s,]) + G$tau.mu[s]
eps.s <- Y.s - mu.s
kappa.hat <- sum(G$theta.kappa * G$X[s,])
xi.s <- sum(G$theta.xi * G$X[s,]) + G$tau.xi[s]

tau.hat.s <- sum(-C.inv[s,-s]/C.inv[s,s] * tau[-s])
varsigma.s <- 1/C.inv[s,s] ##precision matrix stuff



### Testing:
g.s <- g.eta.prime(tau[s], tau.hat.s, varsigma.s, xi.s, log(kappa.hat), eps.s)
gg.s <- g.eta.double.prime(tau[s], tau.hat.s, varsigma.s, xi.s, log(kappa.hat), eps.s)

g.s.arma <- g_eta_prime_new(tau[s], tau.hat.s, varsigma.s, xi.s, log(kappa.hat), eps.s)
gg.s.arma <- g_eta_double_prime_new(tau[s], tau.hat.s, varsigma.s, xi.s, log(kappa.hat), eps.s)


g.s
g.s.arma

gg.s
gg.s.arma



### Testing:
g.s <- g.eta.prime(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)
gg.s <- g.eta.double.prime(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)

g.s.arma <- g_eta_prime_new(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)
gg.s.arma <- g_eta_double_prime_new(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)


g.s
g.s.arma

gg.s
gg.s.arma




ptm=proc.time()
for (i in 1:(5*10^5))
{
  g.s <- g.eta.prime(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)
}
proc.time()-ptm
#user  system elapsed 
#5.212   0.016   5.229

ptm=proc.time()
for (i in 1:(5*10^5))
{
  g.s.arma <- g_eta_prime_new(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)
}
proc.time()-ptm
#user  system elapsed 
#4.976   0.000   4.984 


ptm=proc.time()
for (i in 1:(5*10^5))
{
  gg.s <- g.eta.double.prime(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)
}
proc.time()-ptm
#user  system elapsed 
#5.868   0.008   5.894

ptm=proc.time()
for (i in 1:(5*10^5))
{
  gg.s.arma <- g_eta_double_prime_new(tau[s], tau.hat.s, varsigma.s, 0, log(kappa.hat), eps.s)
}
proc.time()-ptm
#user  system elapsed 
#5.672   0.008   5.696 

### 5 % speedup



