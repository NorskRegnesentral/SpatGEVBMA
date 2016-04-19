rm(list = ls())

library(SpatialGEVBMA)

library(ncdf4io)

setwd("~/NR/SpatGEV")

l.all <- list.files("./inputs/nc_files_used/",pattern = "*.nc")

var.name <- "altCoord"
proj.name="UTM_Zone_33"
proj.att="projection"
topdown=FALSE
xdim.name=NULL
ydim.name=NULL
t=NULL

a <- list()
nm <- NULL
for(i in 1:length(l.all)){
    a[[i]] <- nc4.matrix(paste("./inputs/nc_files_used/",l.all[i],sep=""))  # This file gives x-coordinates on the y-axis and y-coordinates on the x-axis
    nm[i] <- strsplit(l.all[i],".",fixed=TRUE)[[1]][1]
    #X11()
    #image.plot(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i])
    print(i)
  }


# Using the locations from the first file for all of them, assuming they are identical 
# -- which they all are

nx <- length(a[[1]]$x)
ny <- length(a[[1]]$y)

allX <- rep(a[[1]]$x,times=ny)
allY <- rep(a[[1]]$y,each=nx)

allZ <- NULL
b <- a
for (i in 1:length(l.all)){
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

# Saving the grid data here
saveRDS(gridData,file="./inputs/gridData.rds")

# Saving also the data on the original list format
gridDataList <- b
names(gridDataList) <- nm
saveRDS(gridDataList,file="./inputs/gridDataList.rds")
