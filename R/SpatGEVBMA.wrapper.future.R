SpatGEVBMA.wrapper.prediction <- function(mcmc.res, #results file from .inference function.
                                         covariates.folder, # Path to folder with covariate files in netcdf-format (see above) 
                                         output.path = getwd(),  # Path to the where the result folder should be stored
                                         output.folder.name = "SpatGEVBMA.predictions",  # Name of result folder
                                         return.period = c(2,5,10,20,25,50,100,200),  # Return period to impute results for (single number or a vector of numbers)
                                         post.quantiles = c(0.025,0.5,0.975),  # Vector of quantiles for which the posterior should be evaluated
                                         show.uncertainty = TRUE,  # Logical indicating whether an IQR uncertainty plot should also be provided
                                         coordinate.type = "XY", # Character indicating the type/name of coordinate system being used, either "XY" or "LatLon" (see above)
                                         transform.output = NULL, # Character specifying whether and how the output should be transformed. NULL corresponds to no transformation. "UTM_QQ_to_LatLon" transforms from UTM QQ (insert number) to LatLon
                                         burn.in = 10000, # The length of the initial burn-in period which is removed
                                         cores = 1, # The number of cores on the computer used for the imputation. Using detectCores()-1 is good for running on a laptop.
                                         annualMax.name = NULL, # Name of annualMax data used in output plots and netcdf files. If NULL, then the name of the specified sheet is used.
                                         create.tempfiles = FALSE, # Logical indicating whether temporary files should be saved in a Temp folder to perform debugging and check intermediate variables/results if the function crashes
                                         keep.temp.files = FALSE, # Logical indicating whether the temporary files (if written) should be kept or deleted on function completion
                                         save.all.output = TRUE, # Logical indicating whether all R objects should be save to file upon function completion. Allocates approx 2.5 Gb for all of Norway.
                                         testing = FALSE, # Variable indicating whether the run is a test or not. FALSE indicates no testing, a positive number indicates the number of locations being imputed
                                         seed = 123, # The seed used in the mcmc computations
                                         fixed.xi = NULL,  # Where we want the shape parameter fixed
                                         xi.constrain = c(-Inf,Inf))
{
  if (show.uncertainty)
  {
    all.post.quantiles <- c(post.quantiles,c(0.25,0.75))   # The use of sort and unique here messes up things below, so avoid using it.
  }
  
  initial.ls <- ls()  # To be used to subtract globally specified variables when saving intermediate variables
  input.list <- names(formals(SpatGEVBMA.wrapper)) # Want to keep the input variables
  
  
  ## Reading in grid covariates
  cov.files <- list.files(covariates.folder,pattern = "*.nc")
  cov.files.path <- list.files(covariates.folder,pattern = "*.nc",full.names=TRUE)
  
  a <- list()
  nm <- NULL
  for(i in 1:length(cov.files))
  {
    a0 <-   nc_open(cov.files.path[i]) 
    a[[i]] <- list()
    if (coordinate.type=="XY")
    {
      a[[i]]$x <- a0$dim$X$vals
      a[[i]]$y <- a0$dim$Y$vals
    }
    
    if (coordinate.type=="LatLon")
    {
      a[[i]]$x <- a0$dim$Lon$vals
      a[[i]]$y <- a0$dim$Lat$vals
    }
    ## Sorting the variables such that the appear in increasing coordinate order (fields default)
    order.y <- length(a[[i]]$y):1
    a[[i]]$y <- a[[i]]$y[order.y]
    
    a[[i]]$z <- ncvar_get(a0)[,order.y]
    nc_close(a0)
    nm[i] <- strsplit(cov.files[i],".",fixed=TRUE)[[1]][1]
    
    cat(paste("Finished reading ",i," of ",length(cov.files)," covariate files.\n",sep=""))
  }
  
  indX <- a[[1]]$x
  indY <- a[[1]]$y
  
  nx <- length(indX)
  ny <- length(indY)
  
  allX <- rep(indX,times=ny)
  allY <- rep(indY,each=nx)
  
  allZ <- NULL
  b <- a
  for (i in 1:length(cov.files))
  {
    z.vec <- c(a[[i]]$z)
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
  
  gridDataList <- b
  names(gridDataList) <- nm
  
  ## Deleting ununsed large objects
  rm(a,b)
  
  cat("\nCheckpoint 1: Finished structuring of covariate grid.\n")
  
  
  
  ## Checkpoint 3
  R=mcmc.res
  ## Removing burn-in
  R$THETA <- R$THETA[-(1:burn.in),,]
  R$TAU <- R$TAU[-(1:burn.in),,]
  R$ALPHA <- R$ALPHA[-(1:burn.in),]
  R$M <- R$M[-(1:burn.in),,]
  R$LAMBDA <- R$LAMBDA[-(1:burn.in),]
  R$ACCEPT.TAU <- R$ACCEPT.TAU[-(1:burn.in),,]
  #---------------------------------------------#
  
  
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
  
  l_all <- mclapply(1:N, "imputation.func", mc.cores = cores, mc.silent=FALSE,
                    cov.map=cov.map,S.map=S.map,R=R,sigma.22.inv=sigma.22.inv,
                    sigma.22.inv.tau=sigma.22.inv.tau,return.period=return.period,
                    all.post.quantiles=all.post.quantiles,N=N,
                    xi.constrain = xi.constrain)
  l = list()
  l_param = list()
  for(i in 1:length(l_all))
  {
    l[[i]] = l_all[[i]]$Q
    l_param[[i]] = l_all[[i]]$P_Q
  }
  
  Z.p <- array(data = unlist(l),dim = c(length(return.period),length(all.post.quantiles),N))
  Param.maps <- array(data = unlist(l_param),dim = c(3,length(all.post.quantiles),N))
  save(Z.p, Param.maps, file = paste0(output.folder,"/imputation.RData"))
  
  ## START HERE
  
  if (testing)
  {
    N <- dim(cov.map)[1]
    Z.temp <- Z.p
    Z.p <- array(data = NA, dim = c(length(return.period),length(all.post.quantiles), N))
    w.i <- c(1:N0,sample(1:N0,N-N0,replace=TRUE))
    Z.p <- Z.temp[,,w.i,drop=FALSE]
  }
  
  ## Saving intermediate values to easily continue from here if bugs occurs
  if (create.tempfiles){
    current.ls <- ls()
    keep.var <- unique(c(current.ls[!(current.ls %in% initial.ls)],input.list))
    save(list=keep.var,file=file.path(output.temp.folder,"temp_checkpoint_4.RData")) 
    # Here we could delete variables which are not to be used below, to save RAM
  }
  cat("\nCheckpoint 4: Finished mapping posterior to grid.\n")
  
  
  # Checkpoint 4
  
  ## Transforming coordinate system output if applicable
  
  ## Need one of the input files to specify parameters in the ncdf output file
  ncnc <- nc_open(cov.files.path[1])
  
  
  if(!is.null(transform.output)){
    UTM.zone <- as.numeric(substr(transform.output,start=5,stop=6))
    
    # Assume coordinate.type=="XY"
    rangeX <- range(indX)
    rangeY <- range(indY)
    
    allXYMat <- data.frame(X=allX,Y=allY)
    
    #Set projection and zone
    attr(allXYMat, "projection") <- "UTM"
    attr(allXYMat, "zone") <- UTM.zone
    
    #Compute LL coordinates
    allLLMat <- as.matrix(round(convUL(allXYMat, km=FALSE), digits=4))
    
    indLon <- seq(from = min(allLLMat[,1]),to = max(allLLMat[,1]),length.out = nx)
    indLat <- seq(from = min(allLLMat[,2]),to = max(allLLMat[,2]),length.out = ny)
    
    allLon <- rep(indLon,each=ny)
    allLat <- rep(indLat,times=nx)
    
    allLL <- cbind(X=allLon,Y=allLat)
    
    # Set projection and zone
    attr(allLL, "projection") <- "LL"
    attr(allLL, "zone") <- UTM.zone
    
    #Compute UTM coordinates
    XYGrid <- as.matrix(round(convUL(allLL, km=FALSE), digits=0))
    
    # Transform the S matrix to LatLon as well
    #S.new <- cbind(X=S[,1],Y=S[,2])
    
    #attr(S.new, "projection") <- "UTM"
    #attr(S.new, "zone") <- UTM.zone
    
    #Compute LL coordinates
    #S <- as.matrix(round(convUL(S.new, km=FALSE), digits=4))
    
    S=R0$S
  }
  
  ## Writing final results to netcdf-file and image plots
  
  ## Define the dimensions
  if (coordinate.type=="XY")
  {
    if (is.null(transform.output))
    {
      
      output.x <- indX
      output.y <- indY
      
      x.ncdf <- ncdim_def( "X", "meters", output.x)
      y.ncdf <- ncdim_def( "Y", "meters", output.y[ny:1])
      
      dim.list <- list(X=x.ncdf,Y=y.ncdf)
    } else {
      
      original.image <- list(x=indX,y=indY)
      
      output.x <- indLon
      output.y <- indLat
      
      x.ncdf <- ncdim_def( "Lon", "degrees_E", output.x)
      y.ncdf <- ncdim_def( "Lat", "degrees_N", output.y[ny:1])
      
      dim.list <- list(Lon=x.ncdf,Lat=y.ncdf)
      
    }
    
  }
  
  if (coordinate.type=="LatLon")
  {
    output.x <- indX
    output.y <- indY
    
    x.ncdf <- ncdim_def( "Lon", "degrees_E", output.x)
    y.ncdf <- ncdim_def( "Lat", "degrees_N", output.y[ny:1])
    dim.list <- list(Lon=x.ncdf,Lat=y.ncdf)
  }
  
  
  
  ##  Print Return Period Maps 
  for(j in 1:length(return.period))
  {
    ## Just general definitions
    shortName <- paste("quant_",gsub(".","_",post.quantiles,fixed=TRUE),sep="")  # The gsub thing replaces the dot with a underscore
    longName <- paste(post.quantiles," quantile of the marginal ",
                      "posterior distribution for the maximum precipition over ",
                      return.period[j]," years based on data: ",
                      annualMax.name,".",sep="")
    w_median = which(post.quantiles == 0.5)
    if(length(w_median) > 0)
    {
      longName[w_median] <- paste("Median of the marginal posterior ",
                                  "distribution for the ",return.period[j],
                                  " return level for precipitation based on data: ",
                                  annualMax.name,".",sep="")
    }
    
    IQRLongName <- paste("Interquartile range uncertainty measure: Difference" ,
                         "between 0.75-quantile and 0.25-quantile for ",
                         "the maximum precipitaion over ",
                         return.period[j], " years based on data: ",
                         annualMax.name,".",sep="")
    
    filename.nc <- paste0(output.folder,"/posterior.grid_return_",return.period[j],".nc",sep="")
    filename.pdf <- paste0(output.folder,paste("/posterior.return.level.",return.period[j],"grid.pdf",sep=""))
    
    main.quantile = paste("Posterior ", post.quantiles, "-quantile \n ", return.period[j]," year return value with ", annualMax.name," data",sep="")
    main.iqr = paste("Interquartile range uncertainty plot \n ", return.period[j]," year return value with ", annualMax.name," data",sep="")
    output.name <- paste(filename.nc,"_return_",return.period[j],".nc",sep="")
    
    print_map(Q = Z.p[j,,],
              shortName = shortName,
              longName = longName,
              post.quantiles = post.quantiles,
              IQRLongName = IQRLongName,
              show.uncertainty = TRUE,
              ww.na = ww.na,
              n = n,
              output.x = output.x,
              output.y = output.y,
              output.name = output.name,
              filename.nc = filename.nc,
              filename.pdf = filename.pdf,
              main.quantile = main.quantile,
              main.iqr = main.iqr,
              nx = nx,
              ny = ny,
              dim.list = dim.list,
              all.post.quantiles = all.post.quantiles,
              transform.output = transform.output,
              original.image = original.image,
              XYGrid = XYGrid,
              coordinate.type = coordinate.type,
              S = S)
    
  }
  
  ## Print out the parameter maps
  nms_param = c("Location","Inverse Scale","Shape")
  for(j in 1:length(nms_param))
  {
    ## Just general definitions
    shortName <- paste("quant_",gsub(".","_",post.quantiles,fixed=TRUE),sep="")  # The gsub thing replaces the dot with a underscore
    longName <- paste(post.quantiles," quantile of the marginal ",
                      "posterior distribution for the ", nms_param[j], " parameter",
                      " based on data:",
                      annualMax.name,".",sep="")
    
    w_median = which(post.quantiles == 0.5)
    if(length(w_median) > 0)
    {
      longName[w_median] <- paste("Median of the marginal",
                                  "posterior distribution for the ", nms_param[j], " parameter ",
                                  "based on data:",
                                  annualMax.name,".",sep="")
    }
    
    IQRLongName <- paste("Interquartile range uncertainty measure: Difference ",
                         "between 0.75-quantile and 0.25-quantile for ",
                         "the ", nms_param[j], " parameter ",
                         "based on data: ",
                         annualMax.name,".",sep="")
    
    filename.nc <- paste0(output.folder,"/posterior.grid_param_", nms_param[j],".nc",sep="")
    filename.pdf <- paste0(output.folder,"/posterior.param.", nms_param[j],".grid.pdf",sep="")
    
    main.quantile = paste("Posterior ", post.quantiles, "-quantile \n ", "For ", nms_param[j], " parameter with ", annualMax.name," data",sep="")
    main.iqr = paste("Interquartile range uncertainty plot \n ", "For ", nms_param[j], " parameter with ", annualMax.name," data",sep="")
    
    print_map(Q = Param.maps[j,,],
              shortName = shortName,
              longName = longName,
              post.quantiles = post.quantiles,
              IQRLongName = IQRLongName,
              show.uncertainty = TRUE,
              ww.na = ww.na,
              n = n,
              output.x = output.x,
              output.y = output.y,
              output.name = output.name,
              filename.nc = filename.nc,
              filename.pdf = filename.pdf,
              main.quantile = main.quantile,
              main.iqr = main.iqr,
              nx = nx,
              ny = ny,
              dim.list = dim.list,
              all.post.quantiles = all.post.quantiles,
              transform.output = transform.output,
              original.image = original.image,
              XYGrid = XYGrid,
              coordinate.type = coordinate.type,
              S = S)
  }
  
  # Finally save all R objects if desired
  if (save.all.output)
  {
    current.ls <- ls()
    keep.var <- unique(c(current.ls[!(current.ls %in% initial.ls)],input.list))
    save(list=keep.var,file=file.path(output.folder,"final_output.RData")) 
  }
  # Delete temporary files (?)
  if (!keep.temp.files)
  {
    unlink(output.temp.folder,recursive = TRUE)
  }
  
  cat("\nFunction run complete!\n")
  # Function completed!
  
}



