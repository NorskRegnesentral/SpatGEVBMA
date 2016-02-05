
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
Y.s=G$Y.list[[s]]
mu.s <- sum(G$theta.mu * G$X[s,]) + G$tau.mu[s]
kappa.s <- sum(G$theta.kappa * G$X[s,]) + G$tau.kappa[s]
xi.curr <- sum(G$theta.xi * G$X[s,]) + G$tau.xi[s]

### Testing


gg=gev.like(Y.s,mu.s,kappa.s,xi.curr)
gg.arma=as.vector(gev_like_new(Y.s,mu.s,kappa.s,xi.curr))

gg-gg.arma


xi.curr=0
gg=gev.like(Y.s,mu.s,kappa.s,xi.curr)
gg.arma=as.vector(gev_like_new(Y.s,mu.s,kappa.s,xi.curr))

gg-gg.arma


## Identical except for the form of the putput (vector vs. column vector)


### Timing:


ptm=proc.time()
for (i in 1:(5*10^5))
{
  gg=gev.like(Y.s,mu.s,kappa.s,xi.curr)
}
proc.time()-ptm
#user  system elapsed 
#4.540   0.020   5.705 

ptm=proc.time()
for (i in 1:(5*10^5))
{
  gg.arma=gev_like_new(Y.s,mu.s,kappa.s,xi.curr)
}
proc.time()-ptm
#user  system elapsed 
#5.528   0.000   6.279 

## Actually a slowdown  of using C++ here.

