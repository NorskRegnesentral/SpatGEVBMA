rm(list = ls())

library(SpatialGEVBMA)

library(ncdf4io)

setwd("~/NR/SpatGEV")

l.all <- list.files("./inputs/til_Alex/",pattern = "*.nc")

file <- paste("./inputs/til_Alex/",l.all[1],sep="")
var.name <- "altCoord"
proj.name="UTM_Zone_33"
proj.att="projection"
topdown=FALSE
xdim.name=NULL
ydim.name=NULL
t=NULL

a <- list()
nm <- NULL
for(i in 1:length(l.all))
  {
    a[[i]] <- nc4.matrix(paste("./inputs/til_Alex/",l.all[i],sep=""))
    nm[i] <- strsplit(l.all[i],".",fixed=TRUE)[[1]][1]
    X11()
    image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i])
  }



##combine and save as cov.RData
