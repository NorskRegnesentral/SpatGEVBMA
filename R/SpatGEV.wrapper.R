
#### Basic manual and input data file format requirements -- to be put in help file ####

# This is a wrapper for running the spatial.gev.bma function in the SpatGEV package
# and imputing the results on a spatial grid with covariate values, based on a
# netcdf-file for the covariate grid, a spreadsheet with responses at observation stations
# and their spatial location in a text-file. The output is a netcdf-file with user specified
# posterior quantiles for the user specified return period, in addition to spatial maps of 
# these and an optional InterQuartile Range (IQR) uncertainty map. In addtion the input variables
# are written to file.


# The grid covariates values are placed in netcdf-files in the folder 'covariate.folder'.
# The spatial coordinate system can be whatever you like as long as it is decimal spaced 
# (not minutes and seconds) and the same coordinate system is used for the location of 
# the stations (see below). The coordinates must either be specified by captial 'X' and 'Y' 
# (set 'coordinate.type="XY"') or 'Lat' and 'Lon', with capital Ls (set 'coordinate.type="LatLon"').
# All netcdf-files must give covariate values at the same spatial locations (there is currently no check for this)

# The response data for stations are gathered in a spreadsheet located at 'station.annualMax.file',
# at spreadsheet named or numbered 'station.annualMax.sheet'. The first column of this spreadsheet contains
# the observation year, while first row contains the station number (starting from column 2)

# The spatial locations for the station are given in a seperate table formatted txt-file with the first row 
# containing the column names. A column named 'Stnr' must be present to denote the station number corresponding to those
# in the spreadsheet. If 'coordinate.type="XY"', then there must be columns 'X' and 'Y' specifying
# their spatial location in the same coordinate system as used in the netcdf covariate files. 
# If 'coordinate.type="LatLon"', there must be columns 'Lat' and 'Lon' corresponding to the 
# format in the netcdf-file.

# All stations needs to be within the area covered by the grid of covariate values

# Results are imputed only in the locations where all covariates are available


### TO DO ###

## Bugs/quality control
# Ensure that the posterior is correctly written to the nc-file
# Test procedure on data with LatLon coordinates: A) No bugs? B) Does it gives the same results?

## Features
# Allow user specified prior distribution -- currently hard-coded to the default.
# Check for colinearity in the covariates and throw out variables if this 
  # happens -- the current function gives an error on this event and the solution is to 
  # throw out the netcdf file messing it all up from the covariate folder.



SpatGEV.wrapper <- function(covariates.folder, # Path to folder with covariate files in netcdf-format (see above) 
                            station.annualMax.file, # File name of spreadsheet annualMax file (see above)
                            station.annualMax.sheet = 1, # The sheet name or index containing the station annualMax to be read (exactly 1 number)
                            station.locations.file, # File name of table formatted textfile including the spatial locations of the stations 
                            output.path = getwd(),  # Path to the where the result folder should be stored
                            output.folder.name = "SpatGEV.res",  # Name of result folder
                            return.period = 20,  # Return period to impute results for (exactly 1 number)
                            post.quantiles = c(0.025,0.5,0.975),  # Vector of quantiles for which the posterior should be evaluated
                            show.uncertainty = TRUE,  # Logical indicating whether an IQR uncertainty plot should also be provided
                            coordinate.type = "XY", # Character indicating the type/name of coordinate system being used, either "XY" or "LatLon" (see above)
                            table.format = "html", # Character indicating the format for the covariate effect summary tables. Either "html" or "latex".
                            mcmc.reps = 10^5, # Number of MCMC runs for fitting the model with the station data. Should typically be at least be 10^5
                            burn.in = round(mcmc.reps*0.2), # The length of the initial burn-in period which is removed
                            cores = 1, # The number of cores on the computer used for the imputation. Using detectCores()-1 is good for running on a laptop.
                            annualMax.name = NULL, # Name of annualMax data used in output plots and netcdf files. If NULL, then the name of the specified sheet is used.
                            create.tempfiles = FALSE, # Logical indicating whether temporary files should be saved in a Temp folder to perform debugging and check intermediate variables/results if the function crashes
                            keep.temp.files = FALSE, # Logical indicating whether the temporary files (if written) should be kept or deleted on function completion
                            save.all.output = TRUE, # Logical indicating whether all R objects should be save to file upon function completion. Allocates approx 2.5 Gb for all of Norway.
                            testing = FALSE) # Variable indicating whether the run is a test or not. FALSE indicates no testing, a positive number indicates the number of locations being imputed

  ## Various initial fixing
  {
  # Bookeeping for storing intermediate results
  initial.ls <- ls()  # To be used to subtract globally specified variables when saving intermediate variables
  input.list <- names(formals(SpatGEV.wrapper)) # Want to keep the input variables
    
  output.folder <- file.path(output.path,output.folder.name)
  output.temp.folder <- file.path(output.path,output.folder.name,"Temp")

  if (show.uncertainty)
    {
      all.post.quantiles <- sort(unique(c(post.quantiles,c(0.25,0.75))))
    }
  
  ## Initial handling of directories

  # Checks if output directory exists, if not it creates it  
  dirs <- list.dirs(output.path,full.name=FALSE,recursive=FALSE)
  if (!(output.folder.name %in% dirs))
    {
      dir.create(output.folder)
    }

  # The same with the Temp folder 
  dirs0 <- list.dirs(output.folder,full.name=FALSE,recursive=FALSE)
  if (!("Temp" %in% dirs0))
    {
      dir.create(output.temp.folder)
    }

  # Saving the input variables
  save(list=input.list,file=file.path(output.folder,"input_var.RData"))
  
  if (create.tempfiles)
    {
      ## Saving intermediate values to easily continue from here if bugs occurs
      current.ls <- ls()
      keep.var <- unique(c(current.ls[!(current.ls %in% initial.ls)],input.list))
      save(list=keep.var,file=file.path(output.temp.folder,"temp_checkpoint_0.RData")) 
      ## Here we could delete variables which are not to be used below, to save RAM
    }
  
  cat("\nCheckpoint 0: Finished initialization.\n")
  
  # Checkpoint 0
    
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
          a[[i]]$x <- a0$dim$Lat$vals
          a[[i]]$y <- a0$dim$Lon$vals
        }
      ## Sorting the variables such that the appear in increasing coordinate order.
      order.x <- order(a[[i]]$x)
      order.y <- order(a[[i]]$y)
      a[[i]]$x <- a[[i]]$x[order.x]
      a[[i]]$y <- a[[i]]$y[order.y]
      
      a[[i]]$z <- ncvar_get(a0)[order.x,order.y]
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
  
  # Saving the grid data here
  #saveRDS(gridData,file=file.path(output.temp.folder,"gridData.rds"))
  
  # Saving also the data on the original list format
  gridDataList <- b
  names(gridDataList) <- nm
  
  ## Deleting ununsed large objects
  rm(a,b)
  
  if (create.tempfiles)
    {
      ## Saving intermediate values to easily continue from here if bugs occurs
      current.ls <- ls()
      keep.var <- unique(c(current.ls[!(current.ls %in% initial.ls)],input.list))
      save(list=keep.var, file=file.path(output.temp.folder,"temp_checkpoint_1.RData")) 
      ## Here we could delete variables which are not to be used below, to save RAM
    }

  cat("\nCheckpoint 1: Finished structuring of covariate grid.\n")
  
  # Checkpoint 1
  
  ## Reading in station data and extracting corresponding covariates
  fileYData <- loadWorkbook(station.annualMax.file)
  allYData <- suppressWarnings(readWorksheet(fileYData, sheet=station.annualMax.sheet))  # Ignore warnings
  
  SData <- read.table(station.locations.file,header=TRUE)
  if (coordinate.type=="XY")
    {
      x.S <- SData$X
      y.S <- SData$Y
    }
  if (coordinate.type=="LatLon")
    {
      x.S <- SData$Lat
      y.S <- SData$Lon
    }
  SData.use <- data.frame(Stnr=SData$Stnr, x=x.S, y=y.S)

  dat <- suppressWarnings(sapply(allYData,as.numeric)[,-1]) # Ignore warnings
  stations <- as.numeric(substring(colnames(dat),first=2)) # Updated further down
  
  # Extracting Y
  Y.list <- list()
  remove.station <- rep(FALSE,dim(dat)[2])
  for (j in 1:dim(dat)[2])
    { ## Assuming the first column of the Ydata file contains the year
      ## Extracting Ys
      Y.list[[j]] <- c(na.omit(dat[,j]))  # Removing NAs and ignoring the observation year
      if (length(Y.list[[j]]) < 10)
        {
          remove.station[j] = TRUE  # Removes the station if there are less than 10 observations in it
          cat(paste("\nStation ",stations[j]," has only ",length(Y.list[[j]]), " observations, and is therefore removed.\n",sep=""))
        }
    }
  Y.list[remove.station] <- NULL   # Removes station and moves the rest up with single brackets
  
  stations <- as.numeric(substring(colnames(dat),first=2))[!remove.station]  # Assuming first row of Ydata file contains the station numbers (starting from column 2)
  nS <- length(stations)
  
  ## Getting the name from the sheet:
  
  if (is.null(annualMax.name))
    {
      if(is.character(station.annualMax.sheet))
        {
          annualMax.name <- station.annualMax.sheet
        } else {
          annualMax.name <- getSheets(fileYData)[station.annualMax.sheet]
        }
    }
  
                                        # Extracting S
  S <- matrix(NA,nrow=nS,ncol=2)  # Matrix with spatial location of the stations
  
  for (j in 1:nS)
    {
      thisStation <- which(SData.use$Stnr==stations[j])
      S[j,1] <- SData.use$x[thisStation]
      S[j,2] <- SData.use$y[thisStation]
    }
  colnames(S) <- c("x","y")

  ## Extracting X
  # Basic function to be used to pick the closest value when interpolation gives NA values
  get.nn <- function(data, labels, query)
    {
      nns <- get.knnx(data, query, k=1)
      labels[nns$nn.index]
    }
  
  nX=length(gridDataList)
  X = matrix(NA,ncol=nX,nrow=nS)
  for (j in 1:(nX))
    {
      X[,j]=interp.surface(obj=gridDataList[[j]],loc=S)
      theseNA <- which(is.na(X[,j]))
      if (length(theseNA)>0)
        {
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
  

  # Updating X, S and stations after station removal
  X <- cbind(1,X)
  
  nX <- dim(X)[2] # Just updating this one...

  # Putting all station data in a list
  StationData <- list()
  StationData$Y.list <- Y.list
  StationData$X <- X
  StationData$S <- S
  
  
  
  ## Go ahead and save the data lists here as individual files with their names corresponding to the 
  # name of the file.
  
  if (create.tempfiles)
    {
      ## Saving intermediate values to easily continue from here if bugs occurs
      current.ls <- ls()
      keep.var <- unique(c(current.ls[!(current.ls %in% initial.ls)],input.list))
      save(list=keep.var, file=file.path(output.temp.folder,"temp_checkpoint_2.RData")) 
      ## Here we could delete variables which are not to be used below, to save RAM
    }
    
  cat("\nCheckpoint 2: Finished reading station data.\n\n")

  ## Checkpoint 2
  
  ## Running SpatGEV
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
  
  prior$mu$beta.0 <- c(mean(unlist(Y.list)),rep(0, p-1))
  ##prior$mu$Omega.0 <- diag(p)##solve(diag(c(10,rep(100,dim(X.all)[2] - 1))))
  ##prior$kappa$Omega.0 <- diag(p)/1e6##solve(diag(c(100,rep(100,dim(X.all)[2] - 1))))
  ##prior$xi$Omega.0 <- diag(p)##solve(diag(c(100,rep(100,dim(X.all)[2] - 1))))
  
  R0 <- spatial.gev.bma(StationData$Y.list, StationData$X, StationData$S, mcmc.reps, prior, print.every = 1e2)
  
  R <- R0  
  
  ## Removing burn-in
  R$THETA <- R$THETA[-(1:burn.in),,]
  R$TAU <- R$TAU[-(1:burn.in),,]
  R$ALPHA <- R$ALPHA[-(1:burn.in),]
  R$M <- R$M[-(1:burn.in),,]
  R$LAMBDA <- R$LAMBDA[-(1:burn.in),]
  R$ACCEPT.TAU <- R$ACCEPT.TAU[-(1:burn.in),,]
  
  ## Create a table with covariate effects and similar to be written as 
  tbl <- gev.process.results(R)
  rownames(tbl$tbl.mu) <- colnames(X)
  rownames(tbl$tbl.kappa) <- colnames(X)
  rownames(tbl$tbl.xi) <- colnames(X)
  
  tbl.names <- names(tbl)
  write.tables <- paste("tbl$",tbl.names,sep="")
  
  if (table.format=="latex")
    {
      table.format.short <- "tex"
    }  else { 
      table.format.short <- table.format 
    }

  for (i in 1:length(write.tables))
    {
      eval(parse(text=paste("xx <- ",write.tables[i],sep="")))
      xtab <- xtable(xx)   # Trust that the defualt ways to select the number of input variables works fine here
      print.xtable(x=xtab,type=table.format,file=file.path(output.folder,paste("summary_",tbl.names[i],".",table.format.short,sep="")))
    }
  
  if (create.tempfiles)
    {
      ## Saving intermediate values to easily continue from here if bugs occurs
      current.ls <- ls()
      keep.var <- unique(c(current.ls[!(current.ls %in% initial.ls)],input.list))
      save(list=keep.var,file=file.path(output.temp.folder,"temp_checkpoint_3.RData")) 
      ## Here we could delete variables which are not to be used below, to save RAM
    }
  
  cat("\nCheckpoint 3: Finished running SpatGEV.\n\n")
  
  ## Checkpoint 3
  
  ## Mapping posterior to grid

  ## Covariates and locations for the complete grid
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
  
  ##### Additional layer simply for testing purposes

  if (testing)
    {
      N0 <- testing
      N <- N0
    }
  
  l <- mclapply(1:N, "imputation.func", mc.cores = cores, mc.silent=FALSE,
                cov.map=cov.map,S.map=S.map,R=R,sigma.22.inv=sigma.22.inv,
                sigma.22.inv.tau=sigma.22.inv.tau,return.period=return.period,
                all.post.quantiles=all.post.quantiles,N=N)
  
  Z.p <- array(data = unlist(l),dim = c(length(all.post.quantiles),length(return.period),N)) ## TEst this
  
  if (testing)
    {
      N <- dim(cov.map)[1]
      Z.temp <- Z.p
      Z.p <- array(data = NA, dim = c(length(all.post.quantiles), length(return.period), N))
      w.i <- sample(1:dim(Z.temp)[3],N,replace=TRUE)
      Z.p <- Z.temp[,,w.i,drop=FALSE]
    }
  
  ## Saving intermediate values to easily continue from here if bugs occurs
  if (create.tempfiles)
    {
      current.ls <- ls()
      keep.var <- unique(c(current.ls[!(current.ls %in% initial.ls)],input.list))
      save(list=keep.var,file=file.path(output.temp.folder,"temp_checkpoint_4.RData")) 
                                        # Here we could delete variables which are not to be used below, to save RAM
    }
  cat("\nCheckpoint 4: Finished mapping posterior to grid.\n")
  
  
  # Checkpoint 4

  ## Writing final results to netcdf-file and image plots
   
  ## Need one of the input files to specify parameters in the ncdf output file
  ncnc <- nc_open(cov.files.path[1])
  ## Define the dimensions
  if (coordinate.type=="XY")
    {
      x.ncdf <- ncdim_def( "X", "meters", ncnc$dim$X$vals)
      y.ncdf <- ncdim_def( "Y", "meters", ncnc$dim$Y$vals)
      ##    t.ncdf <- ncdim_def( "time", units = "days since 1900-01-01 00:00:00", vals=ncnc$dim$time$vals, unlim=TRUE) # No idea what this time variable does, could probably just skip it.
      dim.list <- list(X=x.ncdf,Y=y.ncdf) #,Time=t.ncdf)
    }
  
  if (coordinate.type=="LatLon")
    {
      x.ncdf <- ncdim_def( "Lat", "degrees", ncnc$dim$Lat$vals)
      y.ncdf <- ncdim_def( "Lon", "degrees", ncnc$dim$Lon$vals)
      ##    t.ncdf <- ncdim_def( "time", units = ncnc$dim$time$units, vals=ncnc$dim$time$vals, unlim=TRUE) # No idea what this time variable does, could probably just skip it.
      dim.list <- list(Lat=x.ncdf,Lon=y.ncdf)#,Time=t.ncdf)
    }
  
  ncvar_defList <- list()

  for(j in 1:length(return.period))
    {
      ## Just general definitions
      shortName <- paste("quant_",post.quantiles,sep="")
      longName <- paste(post.quantiles," quantile of the marginal posterior distribution for the maximum precipition over ",return.period[j]," years based on data: ",
                        annualMax.name,".",sep="")
      
      ## Default variables attributes
      output.unit = "millimeter"  
      output.missval <- -999.99
      output.chunksizes <- c(length(ncnc$dim$X$vals)/3,length(ncnc$dim$Y$vals)/3)
      
                                        # Default global attributes  
      output.references <- "Output from the the SpatGEV.wrapper function in the R-package SpatGEVBMA, developed in Dyrrdal, A. V. et al. (2015)"
      output.referencesRpackage <- "https://github.com/NorskRegnesentral/SpatGEVBMA"
      output.referencesPaper <- "Dyrrdal, A. V., Lenkoski, A., Thorarinsdottir, T. L., & Stordal, F. (2015). Bayesian hierarchical modeling of extreme hourly precipitation in Norway. Environmetrics, 26(2), 89-106."
      output.var.version <- "no.version"  
                                        #output.proj4.string <- "+proj=utm+zone=33+ellps=WGS84"  # MJ: This shouldn't really be hardcoded...
      output.prod.date = substr(Sys.time(), 1, 10)
                                        #output.conventions = "CF-1.4"
      output.institution = "Norwegian Computing Center (Norsk Regnesentral)"
      
      for (i in 1:length(post.quantiles))
        {
          if (post.quantiles[i]==0.5)
            {
              longName[i] <- paste("Median of the marginal posterior distribution for the ",return.period[j]," return level for precipitation based on data: ",annualMax.name,".",sep="")
            }
          ncvar_defList[[i]] <- ncvar_def(name=shortName[i],longname=longName[i],
                                          units=output.unit,
                                          dim=dim.list,
                                          missval=output.missval, # How missing values in input data are defined 
                                          chunksizes = output.chunksizes)
        }
      IQRLongName <- paste("Interquartile range uncertainty measure: Difference between 0.75-quantile and 0.25-quantile for the maximum precipitaion over ",
                           return.period[j]," years based on data: ",annualMax.name,".",sep="")
      if (show.uncertainty)
        {
          this.IQR <- length(post.quantiles)+1
          shortName <- c(shortName,"IQR")
          longName <- c(longName,IQRLongName)
          ncvar_defList[[this.IQR]] <- ncvar_def(name=shortName[this.IQR],longname=longName[this.IQR],
                                                 units=output.unit,
                                                 dim=dim.list,
                                                 missval=output.missval, # How missing values in input data are defined 
                                                 chunksizes = output.chunksizes)
        }

      notNA <- which(!(1:n %in% ww.na))
  
      full.Z.p <- matrix(NA,ncol=n,nrow=length(all.post.quantiles))
      full.Z.p[,notNA] <- Z.p[,j,]
  
      filename.nc <- file.path(output.folder,"posterior.grid")
      outputNc <- nc_create(filename=paste(filename.nc,"_return_",return.period[j],".nc",sep=""),vars=ncvar_defList)
      for (r in 1:length(post.quantiles))
        {
          ncvar_put(outputNc,varid=ncvar_defList[[r]],vals=full.Z.p[r,])  ##check
        }
      if (show.uncertainty)
        {
          IQR <- full.Z.p[this.IQR+1,]-full.Z.p[this.IQR,]
          IQR0 <- IQR + 0 # to save a copy which is not transformed below
          ncvar_put(outputNc,varid=ncvar_defList[[this.IQR]],vals=IQR)
        }
  
      ## Putting global attributes to the nc-file:
      ## Default global attributes  
      ncatt_put(outputNc,varid=0, attname = "references", attval = output.references)
      ncatt_put(outputNc,varid=0, attname = "Rpackage", attval = output.referencesRpackage)
      ncatt_put(outputNc,varid=0, attname = "Paper", attval = output.referencesPaper)
      ncatt_put(outputNc,varid=0, attname = "var.version", attval = output.var.version)
      ##ncatt_put(outputNc,varid=0, attname = "proj4.string", attval = output.proj4.string)
      ncatt_put(outputNc,varid=0, attname = "prod.date", attval = output.prod.date)
      ##ncatt_put(outputNc,varid=0, attname = "conventions", attval = output.conventions)
      ncatt_put(outputNc,varid=0, attname = "institution", attval = output.institution)
  
      nc_close(outputNc)

      if (coordinate.type=="XY"){lab.name=c("x-coord","y-coord")}
      if (coordinate.type=="LatLon"){lab.name=c("Lat","Lon")}
      filename.pdf <- file.path(output.folder,paste("posterior.return.level.",return.period[j],"grid.pdf",sep=""))
      pdf(filename.pdf,width=7, height=7)
      for (r in 1:length(post.quantiles))
        {
          retMat <- matrix(full.Z.p[r,], ncol=ny,nrow=nx)  # This should be correct
          image.plot(indX,indY,retMat,main=paste("Posterior ", post.quantiles[r], "-quantile \n ", return.period[j]," year return value with ", annualMax.name," data",sep=""),xlab=lab.name[1],ylab=lab.name[2])
          points(S[,1],S[,2])
        }
      if (show.uncertainty)
        {
          retMat <- matrix(IQR0, ncol=ny,nrow=nx)  # IQR specified above
          image.plot(indX,indY,retMat,main=paste("Interquartile range uncertainty plot \n ", return.period[j]," year return value with ", annualMax.name," data",sep=""),xlab=lab.name[1],ylab=lab.name[2])
          points(S[,1],S[,2])
        }
      dev.off()
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
  cat("\nNote: Ignore NA coercion warnings.\n\n")
  # Function completed!
}







### Additional help functions

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

imputation.func <- function(i,cov.map,S.map,R,sigma.22.inv,sigma.22.inv.tau,return.period,all.post.quantiles,N)
{

  n.return <- length(return.period)
  n.q <- length(all.post.quantiles)
  Q <- matrix(NA,n.return, n.q)

  X.drop <- cov.map[i,]
  S.drop <- S.map[i,,drop=FALSE]
  P <- gev.impute.params(R, X.drop, S.drop, sigma.22.inv, sigma.22.inv.tau)
  
  for(j in 1:length(return.period))
    {
      z <- gev.z.p(1/return.period[j], P[,1], 1/P[,2], P[,3])
      Q[j,] <- quantile(z, all.post.quantiles)
    }
  
  if(i %% 10 == 0)print(paste(round(i/N*100,2), " % of imputation complete.",sep=""))
  return(Q)
}

