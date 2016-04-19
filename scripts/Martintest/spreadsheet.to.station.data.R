## Assumptions

# Response data for stations are gathered in a spreadsheet with different "types" of data at the different sheets
# All stations are the same for all sheets in the spreadsheet
# Spatial locations given in seperate txt-file
# Covariate values can be found by interpolation on a grid with covariate values

## Format YData:
# First row of Ydata file contains the station numbers (starting from column 2)

## Format txt location data:
# Stnr is the name of the station number in the Spatial data file
# X and Y are the names of the X and Y UTM coordinates of the stations


rm(list = ls())

library(SpatialGEVBMA)
library(ncdf4io)
library(XLConnect)
library(fields)
library(FNN)

setwd("~/NR/SpatGEV")

# Input
folder <- "./inputs/station_data"
folder0 <- "./inputs"
filenameYData <- "AM_allDurations.xlsx"
filenameSData <- "metadata_stations_1hour.txt"
filenameXDataList <- "gridDataList.rds"


excludeSheets = 7:10 # Which of the sheets in the spreadsheet of Y-data to exclude.


### Main script
# Reading in various data
fileYData <- loadWorkbook(file.path(folder,filenameYData))
allYData <- readWorksheet(fileYData, sheet=getSheets(fileYData))  # Ignore warnings
SData <- read.table(file.path(folder,filenameSData),header=TRUE)
gridDataList <- readRDS(file.path(folder0,filenameXDataList))


includedSheets = which(!(1:length(allYData) %in% excludeSheets))
allStationData <- list()
dat <- sapply(allYData[[1]],as.numeric)[,-1] # Ignore warnings
stations <- as.numeric(substring(colnames(dat),first=2))  # Assuming first row of Ydata file contains the station numbers (starting from column 2)
nS <- length(stations)

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

X=cbind(1,X)

# Extracting Y
for (i in includedSheets){
  dat <- sapply(allYData[[i]],as.numeric)[,-1] # Ignore warnings
  allStationData[[i]]<- list()
  allStationData[[i]]$Y.list <- list()
  allStationData[[i]]$X <- X
  allStationData[[i]]$S <- S
  for (j in 1:nS){ # Assuming the first column of th Ydata file contains the year
    # Extracting Ys
    allStationData[[i]]$Y.list[[j]] <- c(na.omit(dat[,j]))  # Removing NAs and ignoring the observation year
  }
  
}

names(allStationData) <- names(allYData)[includedSheets]

## Go ahead and save the data lists here as individual files with their names corresponding to the 
# name of the file.

for (i in includedSheets){
  saveRDS(allStationData[[i]], file=paste("./inputs/processed_input/station_data.",names(allStationData)[i],".rds",sep=""))
}



