rm(list = ls())

library(SpatGEVBMA)
library(parallel)

setwd("~/NR/SpatGEV/")

data(norway)

mc.cores <- 16

load("./inputs/cov.RData")  ## This needs be a direct link to netcdf, in some manner

load("./output/gev.output.bma.RData")

get.sigma.22.inv <- function(R, burn=NULL, odens=1e3)
  {
    n.s <- dim(R$S)[1]
    reps <- dim(R$THETA)[1]
    if (is.null(burn)) burn <- round(reps/10)
    I <- round(seq(burn+1, reps, length=odens))
    sigma.22.inv <- array(dim = c(n.s,n.s,3,length(I)))
    D.S <- make.D(R$S, R$S)
    for (i in 1:length(I))
      {
        print(i)
        it <- I[i]
        for(k in 1:3)
          {
            alpha <- R$ALPHA[it, k]
            lambda <- R$LAMBDA[it, k]
            C <- 1/alpha * exp(-D.S/lambda)
            diag(C) <- diag(C) + 1e-05
            sigma.22.inv[,,k,i] <- solve(C)
        }
      }
    return(sigma.22.inv)
  }

get.sigma.22.inv.tau <- function(R, sigma.22.inv, burn=NULL, odens=1e3)
  {
    n.s <- dim(R$S)[1]
    reps <- dim(R$THETA)[1]
    if (is.null(burn)) burn <- round(reps/10)
    I <- round(seq(burn+1, reps, length=odens))
    sigma.22.inv.tau <- array(dim = c(n.s,3,length(I)))
    for(i in 1:length(I))
      {
        print(i)
        it <- I[i]
        for(k in 1:3)
          {
            sigma.22.inv.tau[,k,i] <- sigma.22.inv[,,k,i] %*% R$TAU[it,,k]
          }
      }
    return(sigma.22.inv.tau)
  }

gev.impute.params <- function (R, X.drop, S.drop, sigma.22.inv, sigma.22.inv.tau, burn = NULL, odens=1e3) 

{
    reps <- dim(R$THETA)[1]
    if (is.null(burn)) burn <- round(reps/10)
    I <- round(seq(burn+1, reps, length=odens))
    P <- matrix(0, length(I), 3)
    D.drop <- make.D(S.drop, R$S)
    MU <- R$THETA[I,,1] %*% X.drop
    KAPPA <- R$THETA[I,,2] %*% X.drop
    XI <- R$THETA[I,,3] %*% X.drop
    for (i in 1:length(I))
      {
        it <- I[i]
        alpha <- R$ALPHA[it, 1]
        lambda <- R$LAMBDA[it, 1]
        sigma.11 <- 1/alpha
        sigma.12 <- 1/alpha * exp(-D.drop/lambda)
        tau.hat <- sigma.12 %*% sigma.22.inv.tau[,1,i]
        varsigma <- sigma.11 - sigma.12 %*% sigma.22.inv[,,1,i] %*% t(sigma.12)
        tau.new <- rnorm(1, tau.hat, sd = sqrt(varsigma))
        
        mu.s <- MU[i] + tau.new
        
        alpha <- R$ALPHA[it, 2]
        lambda <- R$LAMBDA[it, 2]
        
        sigma.11 <- 1/alpha
        sigma.12 <- 1/alpha * exp(-D.drop/lambda)
        tau.hat <- sigma.12 %*% sigma.22.inv.tau[,2,i]
        varsigma <- sigma.11 - sigma.12 %*% sigma.22.inv[,,2,i] %*% t(sigma.12)
        tau.new <- rnorm(1, tau.hat, sd = sqrt(varsigma))
        
        kappa.hat <- KAPPA[i]
        kappa.s <- rtnorm(1, kappa.hat + tau.hat, sd = sqrt(varsigma),lower = 0)
        
        alpha <- R$ALPHA[it, 3]
        lambda <- R$LAMBDA[it, 3]
        
        sigma.11 <- 1/alpha
        sigma.12 <- 1/alpha * exp(-D.drop/lambda)
        tau.hat <- sigma.12 %*% sigma.22.inv.tau[,3,i]
        varsigma <- sigma.11 - sigma.12 %*% sigma.22.inv[,,3,i] %*% t(sigma.12)
        tau.new <- rnorm(1, tau.hat, sd = sqrt(varsigma))
        
        xi.s <- XI[i] + tau.new
        
        P[i,] <- c(mu.s, kappa.s, xi.s)
      }
    return(P)
  }


sigma.22.inv <- get.sigma.22.inv(R)
sigma.22.inv.tau <- get.sigma.22.inv.tau(R, sigma.22.inv)

##------------- This fixes stupidity -----------------------
ww.na <- which(apply(is.na(cov),1,"any"))
cov <- cov[-ww.na,]
               
## Can we leave this at "long story"
X.orig <- as.matrix(read.delim("./inputs/mc.dsgn.mat.txt",sep=" "))
mu.orig <- colMeans(X.orig, na.rm=TRUE)
sd.orig <- sqrt(diag(var(X.orig)))
                
cov.map <- cbind(1,cov)

colnames(cov.map) <- c("","lat","lon", "JJAtemp","elev", "distancSea", "MAP", "MSP", "wetDays", "JJAtemp.1")
for(i in 2:dim(cov.map)[2])
  {
    cov.map[,i] <- (cov.map[,i] - mu.orig[i - 1])/sd.orig[i - 1]
  }
S.map <- cov.map[,3:2]
##----------------------------------------------------------

helper <- function(i)
  {
    X.drop <- cov.map[i,]
    S.drop <- S.map[i,,drop=FALSE]
    P <- gev.impute.params(R, X.drop, S.drop, sigma.22.inv, sigma.22.inv.tau)
    z <- gev.z.p(1/20, P[,1], 1/P[,2], P[,3])
    if(i %% 10 == 0)print(paste("Finished",i))
    return(quantile(z, c(.025,.5,.975)))
  }

N <- dim(cov.map)[1]

l <- mclapply(1:N, "helper", mc.cores = mc.cores, mc.silent=FALSE)
Z.p <- matrix(unlist(l),ncol=3,byrow=TRUE)
save(l, file="./output/returns.RData")
