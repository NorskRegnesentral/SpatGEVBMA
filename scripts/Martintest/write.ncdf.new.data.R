rm(list = ls())

library(SpatialGEVBMA)

setwd("~/NR/SpatGEV/")

## Rely on the ncdf4 package in R
library(ncdf4)

quants = c(0.025,0.5,0.975)
subset=FALSE


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

###

### Just need this part for now to extract the subset data below

if (subset){
gridData <- readRDS("./inputs/gridData.rds")
cov <- as.matrix(gridData$covariates)
cov <- cov[1:100000,]
ww.na <- which(apply(is.na(cov),1,"any"))
S.map0 <- as.matrix(gridData$coordinates[1:100000,])
S.map <- S.map0[-ww.na,]
}


#####


# Need one of the input files to specify parameters in the ncdf output file
l.all <- list.files("./inputs/nc_files_used/",pattern = "*.nc")
ncnc <- nc_open(paste("./inputs/nc_files_used/",l.all[1],sep=""))
#cc = ncvar_get(ncnc) # Not needed?

### Just general definitions
# Define the dimensions
x.ncdf <- ncdim_def( "X", "meters", ncnc$dim$X$vals)
y.ncdf <- ncdim_def( "Y", "meters", ncnc$dim$Y$vals)

if (subset){
x.ncdf <- ncdim_def( "X", "meters", unique(S.map0[,1]))
y.ncdf <- ncdim_def( "Y", "meters", unique(S.map0[,2]))
}

t.ncdf <- ncdim_def( "time", units = ncnc$dim$time$units, vals=ncnc$dim$time$vals, unlim=TRUE) # No idea what this time variable does, could probably just skip it.
# I do skip the nb2 variable

ncvar_defList <- list()
shortName <- paste("quant_",quants,sep="")
longName <- paste(quants," quantile of the marginal posterior distribution for the maximum hourly precipition over the year")


for (i in 1:length(quants)){
  if (quants[i]==0.5){
    longName[i] <- "Median of the marginal posterior distribution for the maximum hourly precipition over the year"
  }
  if (subset){
    ncvar_defList[[i]] <- ncvar_def(name=shortName[i],longname=longName[i],
                    units="mm",
                    dim=list(X=x.ncdf,Y=y.ncdf,Time=t.ncdf),
                    missval=-999.99, # How missing values in input data are defined 
                    chunksizes = c(1,1,1))
  } else {
    ncvar_defList[[i]] <- ncvar_def(name=shortName,longname=longName,
                                  units="mm",
                                  dim=list(X=x.ncdf,Y=y.ncdf,Time=t.ncdf),
                                  missval=-999.99, # How missing values in input data are defined 
                                  chunksizes = ncnc$var$precipitation_amount$chunksizes)
  }
}


outputFiles <- list.files("./output/new",pattern = "mapreturns.")
outputFilesFull <- paste("./output/new/",outputFiles,sep="")
fileNames <- getstr(outputFiles,"\\.","\\.")


for (i in 1:length(outputFiles)){
  returns <- readRDS(outputFilesFull[i])
  returnsMat <- matrix(unlist(returns),ncol=length(quants),byrow=TRUE)
  
  filefile <- paste("./output/new/final_returns.",fileNames[i],".nc",sep="")
  outputNc <- nc_create(filename=filefile,vars=ncvar_defList)
  for (j in 1:length(quants)){
    
    if (subset){
    # Just for now:
    valval <- rep(NA,length(unique(S.map0[,1]))*length(unique(S.map0[,2])))
    notNA <- which(!(1:100000 %in% ww.na))
    valval[notNA] <- returnsMat[,j]
    } else{
      valval <- returnsMat[,j]
    }
    ncvar_put(outputNc,varid=ncvar_defList[[j]],vals=valval)
  }
  nc_close(outputNc)
}



### Here I would like to read the nc file, plot it on a map and also include the observations and the returns levels there

library(ncdf4io)

stationData <- readRDS(file="./inputs/processed_input/station_data.AM_10min.rds")

thisNc <- nc_open(filename="./output/new/final_returns.AM_10min.nc")
varNames <- names(thisNc$var)

gridReturns <- list()

for (i in 1:length(varNames)){
  gridReturns[[i]] <- list()
  gridReturns[[i]]$X <- thisNc$dim$X$vals
  gridReturns[[i]]$Y <- thisNc$dim$Y$vals
  gridReturns[[i]]$Z <- ncvar_get(thisNc,varid=varNames[i])
}

for (i in 1:length(varNames)){
  x11()
  image.plot(gridReturns[[i]]$X,gridReturns[[i]]$Y,gridReturns[[i]]$Z,main=paste("AM_10min: ",varNames[i],sep=""))
}



meanprec= rep(NA,length(stationData$Y.list))
for (i in 1:length(stationData$Y.list)){
  meanprec[i] <- mean(stationData$Y.list[[i]])
}

this <- nc4.matrix("./output/new/final_returns5.AM_10min.nc")


points(stationData$S[,1],stationData$S[,2],pch=1:dim(stationData$S)[1])

im=as.image(cov[,6],x=S.map)
image.plot(im)





