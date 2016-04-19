rm(list = ls())

library(SpatialGEVBMA)
library(parallel)

setwd("~/NR/SpatGEV/")

subset=FALSE

#data(norway)

## Just definitions here

getstr = function(mystring, initial.character="_", final.character="_")
{
  # check that all 3 inputs are character variables
  if (!is.character(mystring))
  {
    stop('The parent string must be a character variable.')
  }
  
  if (!is.character(initial.character))
  {
    stop('The initial character must be a character variable.')
  }
  
  
  if (!is.character(final.character))
  {
    stop('The final character must be a character variable.')
  }
  
  add=0
  if(initial.character==final.character){add=1}
  
  # pre-allocate a vector to store the extracted strings
  snippet = rep(0, length(mystring))
  
  for (i in 1:length(mystring))
  {
    # extract the initial position
    initial.position = gregexpr(initial.character, mystring[i])[[1]][1] + 1
    
    # extract the final position
    final.position = gregexpr(final.character, mystring[i])[[1]][1+add] - 1
    
    # extract the substring between the initial and final positions, inclusively
    snippet[i] = substr(mystring[i], initial.position, final.position)
  }
  return(snippet)
}

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

# Number of cores for parallelization
mc.cores <- 20

# Covariates and locations for the complete grid
gridData <- readRDS("./inputs/gridData.rds")
cov <- as.matrix(gridData$covariates)

if (subset){
cov <- cov[1:100000,]
}

ww.na <- which(apply(is.na(cov),1,"any"))
cov <- cov[-ww.na,]
cov.map <- cbind(1,cov)
colnames(cov.map) <- c("",names(gridData$covariates))

S.map <- as.matrix(gridData$coordinates)
if(subset){
  S.map=S.map[1:100000,]
}
S.map <- S.map[-ww.na,]


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

outputFiles <- list.files("./output/new",pattern = "output.")
outputFilesFull <- paste("./output/new/",outputFiles,sep="")
fileNames <- getstr(outputFiles,"\\.","\\.")

for (i in 1:length(outputFiles)){
  load(outputFilesFull[i])
  sigma.22.inv <- get.sigma.22.inv(R)
  sigma.22.inv.tau <- get.sigma.22.inv.tau(R, sigma.22.inv)

  l <- mclapply(1:N, "helper", mc.cores = mc.cores, mc.silent=FALSE)  
  Z.p <- matrix(unlist(l),ncol=3,byrow=TRUE)
  saveRDS(l,file=paste("./output/new/mapreturns.",fileNames[i],".rds",sep=""))
}



