# This is a SpatGEV wrapper taking 
# covariate input through netcdf-files, 
# station returns as a xlsx-file,
# and station locations as a txt-file,
# then running these for the full set of locations specified in the netcdf-files,
# and writing a netcdf-file and some basic figures to a specified folder

library(SpatialGEVBMA)

library(ncdf4io)
library(XLConnect)
library(FNN)
library(parallel)
library(ncdf4)
library(fields)


### TO DO:

## Bugs/quality control
# Check/ask Alex what causes the difference in indX and ncnc$dim$X$vals. What is done in nc4.matrix just seems wrong...
# Find the source of the very large return values we get, there ought to be something wrong there.
# Determine if the posterior is correctly written to the nc-file
# Test the procedure on the norway-data (?)

## Features:
# Specify coordinate system (UTM no or lat/lon)
#X Write basic IQR figure to file
# Document input variables
# Document the form of the input files
# Ensure output is written correctly




SpatGEV.wrapper <- function(covariates.folder, station.returns.file, station.returns.sheet, station.locations.file, output.location, output.folder.name, keep.temp.files,
                            mcmc.reps,cores, post.quantiles, uncertainty.plot, testing,
                            returns.name = NULL){

  # For simplicity setting default values here (to be removed)
  covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used"
  station.returns.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
  station.returns.sheet <- 1 # The sheet name or index containing the station returns to be read
  station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
  output.location <- "~/NR/SpatGEV"
  output.folder.name <- "outputTemp2"
  keep.temp.files <- TRUE
  mcmc.reps <- 5*10^4 # Should at least be 10^4
  cores <- 10 # 20 # The number of cores to be used for paralellization when mapping out the posterior to the grid
  post.quantiles <- c(0.025,0.5,0.975)
  uncertainty.plot <- TRUE
  testing <- FALSE
  returns.name <- NULL # If NULL then the name of specified sheet is used
  
  ### The actual function

  ## Various initial fixing
  {
  output.folder <- file.path(output.location,output.folder.name)
  output.temp.folder <- file.path(output.location,output.folder.name,"Temp")

  if (uncertainty.plot){
    all.post.quantiles <- unique(c(post.quantiles,c(0.25,0.75)))
  }
  
  }
  
  ## Initial handling of directories
  {
  # Checks if output directory exists, if not it creates it  
  dirs <- list.dirs(output.location,full.name=FALSE,recursive=FALSE)
  if (!(output.folder.name %in% dirs)){
    dir.create(output.folder)
  }

  # The same with the Temp folder 
  dirs0 <- list.dirs(output.folder,full.name=FALSE,recursive=FALSE)
  if (!("Temp" %in% dirs0)){
    dir.create(output.temp.folder)
  }
  
  }
  
  ## Help functions
  {
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
        #print(i)
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
        #print(i)
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

  helper <- function(i)
  {
    X.drop <- cov.map[i,]
    S.drop <- S.map[i,,drop=FALSE]
    P <- gev.impute.params(R, X.drop, S.drop, sigma.22.inv, sigma.22.inv.tau)
    z <- gev.z.p(1/20, P[,1], 1/P[,2], P[,3])
    if(i %% 10 == 0)print(paste("Finished",i))
    return(quantile(z, all.post.quantiles))
  }
  }  
    
  ## Reading in grid covariates
  {
  cov.files <- list.files(covariates.folder,pattern = "*.nc")
  cov.files.path <- list.files(covariates.folder,pattern = "*.nc",full.names=TRUE)
  
  a <- list()
  nm <- NULL
  for(i in 1:length(cov.files)){
    a[[i]] <- nc4.matrix(cov.files.path[i])  # This file gives x-coordinates on the y-axis and y-coordinates on the x-axis
    nm[i] <- strsplit(cov.files[i],".",fixed=TRUE)[[1]][1]
    
    cat(paste("Finished reading ",i," of ",length(cov.files)," covariate files.\n",sep=""))
    }
  
  indX <- a[[1]]$x
  indY <- a[[1]]$y
  indZ <- a[[1]]$z

  nx <- length(indX)
  ny <- length(indY)
      
  allX <- rep(indX,times=ny)
  allY <- rep(indY,each=nx)
  
  allZ <- NULL
  b <- a
  for (i in 1:length(cov.files)){
    z.vec <- c(indZ)
    mu.z.vec <- mean(z.vec, na.rm=TRUE)
    sd.z.vec <- sd(z.vec,na.rm=TRUE)  
    stand.z.vec <- (z.vec-mu.z.vec)/sd.z.vec
    allZ <- cbind(allZ,stand.z.vec)
    b[[i]]$z <- matrix(stand.z.vec,ncol=ny)
  }
  
  gridData <- list()
  gridData$coordinates <- data.frame(x=allX,y=allY)
  gridData$covariates <- as.data.frame(allZ)
  colnames(gridData$covariates) <- nm
  gridData$n <- length(allX)
  n <- gridData$n
  
  # Saving the grid data here
  saveRDS(gridData,file=file.path(output.temp.folder,"gridData.rds"))
  
  # Saving also the data on the original list format
  gridDataList <- b
  names(gridDataList) <- nm
  saveRDS(gridDataList,file=file.path(output.temp.folder,"gridDataList.rds"))
  
  cat("\nFinished structuring of covariate grid.\n")
  }
  
  ## Reading in station data
  {
  fileYData <- loadWorkbook(station.returns.file)
  allYData <- readWorksheet(fileYData, sheet=station.returns.sheet)  # Ignore warnings
  
  SData <- read.table(station.locations.file,header=TRUE)
  
  StationData <- list()
  dat <- sapply(allYData,as.numeric)[,-1] # Ignore warnings
  stations <- as.numeric(substring(colnames(dat),first=2))  # Assuming first row of Ydata file contains the station numbers (starting from column 2)
  nS <- length(stations)
  
  # Getting the name from the sheet:
  
  if (is.null(returns.name)){
    if(is.character(station.returns.sheet)){
      returns.name <- station.returns.sheet
    } else {
      returns.name <- getSheets(fileYData)[station.returns.sheet]
    }
  }
  
  # Extracting S
  S <- matrix(NA,nrow=nS,ncol=2)  # Matrix with spatial location of the stations
  
  for (j in 1:nS){
    thisStation <- which(SData$Stnr==stations[j]) # Assuming Stnr is the name of the station number in the Spatial data file
    S[j,1] <- SData$X[thisStation]
    S[j,2] <- SData$Y[thisStation]
  }
  colnames(S) <- c("x","y")

  ## Extracting X
  # Basic function to be used to pick the closetest value when interpolation gives NA values
  get.nn <- function(data, labels, query) {
    nns <- get.knnx(data, query, k=1)
    labels[nns$nn.index]
  }
  
  nX=length(gridDataList) 
  X = matrix(NA,ncol=nX,nrow=nS)
  for (j in 1:(nX)){
    X[,j]=interp.surface(obj=gridDataList[[j]],loc=S)
    theseNA <- which(is.na(X[,j]))
    if (length(theseNA)>0){
      nx <- length(gridDataList[[1]]$x)
      ny <- length(gridDataList[[1]]$y)
      xyMat <- cbind(x=rep(gridDataList[[j]]$x,times=ny),y=rep(gridDataList[[j]]$y,each=nx))
      labs <- c(gridDataList[[j]]$z)
      
      labs <- labs[which(!is.na(labs))]
      xyMat <- xyMat[which(!is.na(labs)),]
      X[theseNA,j] <- get.nn(data=xyMat,labels=labs,query=cbind(x=S[theseNA,1],y=S[theseNA,2]))
    }
  }
  colnames(X) <- names(gridDataList)
  
  X <- cbind(1,X)
  nX <- dim(X)[2] # Just updating this one...
  
  # Extracting Y
  StationData$Y.list <- list()
  StationData$X <- X
  StationData$S <- S
  for (j in 1:nS){ # Assuming the first column of th Ydata file contains the year
    # Extracting Ys
    StationData$Y.list[[j]] <- c(na.omit(dat[,j]))  # Removing NAs and ignoring the observation year
  }
    
  ## Go ahead and save the data lists here as individual files with their names corresponding to the 
  # name of the file.
  
  saveRDS(StationData, file=file.path(output.temp.folder,"stationData.rds"))
  
  cat("\nFinished reading station data.\nNote: Ignore NA coercion warnings.\n")
  }
  
  ## Running SpatGEV
  {
  p <- nX
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
  
  R <- spatial.gev.bma(StationData$Y.list, StationData$X, StationData$S, mcmc.reps, prior, print.every = 1e2)
  tbl <- gev.process.results(R)
  save(R, tbl, file=file.path(output.temp.folder,"mcmc.output.RData"))
  
  cat("\nFinished running SpatGEV.\n")
  }
  
  ## Mapping posterior to grid
  {

  # Covariates and locations for the complete grid
  cov <- as.matrix(gridData$covariates)
  
  ww.na <- which(apply(is.na(cov),1,"any"))
  cov <- cov[-ww.na,]
  cov.map <- cbind(1,cov)
  colnames(cov.map) <- c("",names(gridData$covariates))
  
  S.map <- as.matrix(gridData$coordinates)
  S.map <- S.map[-ww.na,]
  
  N <- dim(cov.map)[1]
  
  sigma.22.inv <- get.sigma.22.inv(R)
  sigma.22.inv.tau <- get.sigma.22.inv.tau(R, sigma.22.inv)
  
  ##### TESTING Simply for testing purposes

  if (testing){
    N0=1000
    N=N0
  }
  l <- mclapply(1:N, "helper", mc.cores = cores, mc.silent=FALSE)  
  Z.p <- matrix(unlist(l),length(all.post.quantiles))
  if (testing){
    N <- dim(cov.map)[1]
    copyMat <- t(matrix(rep(Z.p[,N0],each=N-N0),ncol=length(all.post.quantiles)))
    Z.p <- cbind(Z.p,copyMat)
    }
  saveRDS(Z.p,file=file.path(output.temp.folder,"mapreturns.rds"))
  
    
  cat("\nFinished mapping posterior to grid.\n")
  
  
  }

  ## Writing final results to netcdf-file and image plots
  {
  
  # Need one of the input files to specify parameters in the ncdf output file
  ncnc <- nc_open(cov.files.path[1])


  ### Just general definitions
  # Define the dimensions
  x.ncdf <- ncdim_def( "X", "meters", ncnc$dim$X$vals)
  y.ncdf <- ncdim_def( "Y", "meters", ncnc$dim$Y$vals)
  t.ncdf <- ncdim_def( "time", units = ncnc$dim$time$units, vals=ncnc$dim$time$vals, unlim=TRUE) # No idea what this time variable does, could probably just skip it.

  ncvar_defList <- list()
  shortName <- paste("quant_",post.quantiles,sep="")
  longName <- paste(post.quantiles," quantile of the marginal posterior distribution for the maximum hourly precipition over the year")
  
  
  for (i in 1:length(post.quantiles)){
    if (post.quantiles[i]==0.5){
      longName[i] <- "Median of the marginal posterior distribution for the maximum hourly precipition over the year"
    }
    ncvar_defList[[i]] <- ncvar_def(name=shortName[i],longname=longName[i],
                                      units="mm",
                                      dim=list(X=x.ncdf,Y=y.ncdf,Time=t.ncdf),
                                      missval=-999.99, # How missing values in input data are defined 
                                      chunksizes = ncnc$var$precipitation_amount$chunksizes)
    }
  if (uncertainty.plot){
    this.IQR <- length(post.quantiles)+1
    shortName <- c(shortName,"IQR")
    longName <- c(longName,"Interquartile range uncertainty measure: Difference between 0.75-quantile and 0.25-quantile.")
    ncvar_defList[[this.IQR]] <- ncvar_def(name=shortName[this.IQR],longname=longName[this.IQR],
                                    units="mm",
                                    dim=list(X=x.ncdf,Y=y.ncdf,Time=t.ncdf),
                                    missval=-999.99, # How missing values in input data are defined 
                                    chunksizes = ncnc$var$precipitation_amount$chunksizes)
  }
  
    
  notNA <- which(!(1:n %in% ww.na))
  
  full.Z.p <- matrix(NA,ncol=n,nrow=length(all.post.quantiles))
  full.Z.p[,notNA] <- Z.p
  
  filename.nc <- file.path(output.folder,"posterior.grid.nc")
  outputNc <- nc_create(filename=filename.nc,vars=ncvar_defList)
  for (j in 1:length(post.quantiles)){
    ncvar_put(outputNc,varid=ncvar_defList[[j]],vals=full.Z.p[j,])
  }
  if (uncertainty.plot){
    IQR <- full.Z.p[this.IQR+1,]-full.Z.p[this.IQR,]
    ncvar_put(outputNc,varid=ncvar_defList[[this.IQR]],vals=IQR)
  }
  
  nc_close(outputNc)
  
  filename.pdf <- file.path(output.folder,"posterior.return.grid.pdf")
  pdf(filename.pdf,width=7, height=7)
  for (j in 1:length(post.quantiles)){
    retMat <- matrix(full.Z.p[j,], ncol=ny,nrow=nx)  # This should be correct
    image.plot(indX,indY,retMat,main=paste(returns.name,": Posterior return, ", post.quantiles[j]," quantile",sep=""),xlab="x",ylab="y")
    points(S[,1],S[,2])
  }
  if (uncertainty.plot){
  retMat <- matrix(IQR, ncol=ny,nrow=nx)  # IQR specified above
  image.plot(indX,indY,retMat,main=paste(returns.name,": Interquartile range uncertainty plot",sep=""),xlab="x",ylab="y")
  points(S[,1],S[,2])
  }
  dev.off()
  
  # Should possibly improve this image, adding the possibility to return also an uncertainty plot with the IQR or so.
  
  }
  
}
  
  