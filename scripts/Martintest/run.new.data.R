##  This script replicates the results from Dyrrdal et al.
rm(list = ls())

library(SpatialGEVBMA)

setwd("~/NR/SpatGEV/")

### A helper function


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

### Loading rds files in the inputs/processed_input folder:

stationFilesVec <- list.files("./inputs/processed_input/",pattern = "*.rds")
YName <- getstr(stationFilesVec,'\\.','\\.')

stationFiles <- list()
for (i in 1:length(stationFilesVec)){
  stationFiles[[i]] <- readRDS(paste("./inputs/processed_input/",stationFilesVec[i],sep=""))
}

n.reps <- 1e4

for (i in 1:length(stationFiles)){
  Y <- stationFiles[[i]]$Y
  X <- stationFiles[[i]]$X
  S <- stationFiles[[i]]$S
  p <- dim(X)[2]
  prior <- NULL
  prior$mu$alpha.a <- 2
  prior$mu$alpha.b <- 6
  prior$mu$lambda.a <- 2
  prior$mu$lambda.b <- 2

  prior$kappa$alpha.a <- 2
  prior$kappa$alpha.b <- 2
  prior$kappa$lambda.a <- 1.5
  prior$kappa$lambda.b <- 1.5

  prior$xi$alpha.a <- 2
  prior$xi$alpha.b <- 1
  prior$xi$lambda.a <- 2
  prior$xi$lambda.b <- 1

  prior$mu$beta.0 <- c(8,rep(0, p-1))
  ##prior$mu$Omega.0 <- diag(p)##solve(diag(c(10,rep(100,dim(X.all)[2] - 1))))
  ##prior$kappa$Omega.0 <- diag(p)/1e6##solve(diag(c(100,rep(100,dim(X.all)[2] - 1))))
  ##prior$xi$Omega.0 <- diag(p)##solve(diag(c(100,rep(100,dim(X.all)[2] - 1))))

R <- spatial.gev.bma(Y, X, S, n.reps, prior, print.every = 1e2)
tbl <- gev.process.results(R)
save(R, tbl, file=paste("./output/new/output.",YName[i],".RData",sep=""))
}


