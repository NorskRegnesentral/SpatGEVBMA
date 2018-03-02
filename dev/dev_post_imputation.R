rm(list = ls())

library(SpatGEVBMA)

setwd("~/Desktop/toAlex/SpatGEV.res.1440min/")

output.temp.folder  = "./Temp/"

load(paste0("input_var.RData"))
load(paste0(output.temp.folder, "/temp_checkpoint_0.RData"))
load(paste0(output.temp.folder, "/temp_checkpoint_1.RData"))
load(paste0(output.temp.folder, "/temp_checkpoint_2.RData"))
load(paste0(output.temp.folder, "/temp_checkpoint_3.RData"))
load("./mcmc.RData")
R = R0
load("./imputation.RData")

  
  ## Checkpoint 3
  
  ## Mapping posterior to grid

  ## Covariates and locations for the complete grid
  cov <- as.matrix(gridData$covariates)
  
  ww.na <- which(apply(is.na(cov),1,"any"))
  cov <- cov[-ww.na,]
  cov.map <- cbind(1,cov)
  colnames(cov.map) <- c("",names(gridData$covariates))

  if(coordinate.type == "XY")
  {
    S.map <- as.matrix(gridData$coordinates) / 1e4 ## See above where StationData$S is formed
  }else{
    S.map <- as.matrix(gridData$coordinates)
  }
  S.map <- S.map[-ww.na,]
  
  N <- dim(cov.map)[1]
  
  sigma.22.inv <- get.sigma.22.inv(R)
  sigma.22.inv.tau <- get.sigma.22.inv.tau(R, sigma.22.inv)
  
  ##### Additional layer simply for testing purposes

  if (testing)
  {
    N0 <- testing
    N <- N0
  }
  
  RNGkind("L'Ecuyer-CMRG") # In order to get the same 
  



