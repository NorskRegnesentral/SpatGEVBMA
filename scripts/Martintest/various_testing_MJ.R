#### Reading and writing netcdf files:


### Testing to read one of the covariate files:
library(SpatialGEVBMA)

setwd("~/NR/SpatGEV/")

l.all <- list.files("./inputs/nc_files_used/",pattern = "*.nc")

i=1
aa=nc4.matrix(paste("./inputs/nc_files_used/",l.all[i],sep=""))
bb=nc4in(file=paste("./inputs/nc_files_used/",l.all[i],sep=""))  # Not working

library(ncdf4)

ncnc <- nc_open(paste("./inputs/nc_files_used/",l.all[i],sep=""))
cc = ncvar_get(ncnc)


##### Trying to write some data to netcdf file with teh ncdf4 package:

# dimension definitions:

# Define some straightforward dimensions
x.ncdf <- ncdim_def( "X", "meters", ncnc$dim$X$vals)
y.ncdf <- ncdim_def( "Y", "meters", ncnc$dim$Y$vals)
t.ncdf <- ncdim_def( "time", units = ncnc$dim$time$units, vals=ncnc$dim$time$vals, unlim=TRUE) # No idea what this time variable does, could probably just skip it.
# I do skip the nb2 variable

# Should include the actual 
lowerQ <- ncvar_def(name="l_post_quantile",longname="0.025 (lower) quantile of the marginal posterior distribution",
                units="FILL IN UNIT HERE",
                dim=list(X=x.ncdf,Y=y.ncdf,Time=t.ncdf),
                missval=-999.99, # What 
                chunksizes = ncnc$var$precipitation_amount$chunksizes)
med <- ncvar_def(name="med_post_quantile",longname="Median of the marginal posterior distribution",
                  units="FILL IN UNIT HERE",
                  dim=list(X=x.ncdf,Y=y.ncdf,Time=t.ncdf),
                  missval=-999.99, # What 
                  chunksizes = ncnc$var$precipitation_amount$chunksizes)
upperQ <- ncvar_def(name="u_post_quantile",longname="0.975 (upper) quantile of the marginal posterior distribution",
                  units="FILL IN UNIT HERE",
                  dim=list(X=x.ncdf,Y=y.ncdf,Time=t.ncdf),
                  missval=-999.99, # What 
                  chunksizes = ncnc$var$precipitation_amount$chunksizes)

outputNc <- nc_create(filename="TESTING2.nc",vars=list(lowerQ,med,upperQ))

this=ncvar_get(outputNc,varid="u_post_quantile")

# Now putting data into the ncdf-file

ncvar_put(outputNc,varid=lowerQ,vals=1:(1195*1550))
ncvar_put(outputNc,varid=med,vals=1000*(1:(1195*1550)))
ncvar_put(outputNc,varid=upperQ,vals=rep(2,1195*1550))

this=ncvar_get(outputNc,varid="l_post_quantile")
this=ncvar_get(outputNc,varid="med_post_quantile")


nc_close(outputNc)








aa=NULL
for (i in 1:21){
  aa=cbind(aa,a[[i]]$x)
}

i=1
i=2

x11()
a[[i]] <- nc4.matrix(paste("./inputs/til_Alex/",l.all[i],sep=""))
nm[i] <- strsplit(l.all[i],".",fixed=TRUE)[[1]][1]
image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i])
points(S)
points(S[which(is.na(X[,14])),],col="green")

x11()
image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i],xlim=c(-10^5,0),ylim=c(6400000,6700000))
points(S)
points(S[which(is.na(X[,14])),],col="green")



i=14

x11()
a[[i]] <- nc4.matrix(paste("./inputs/til_Alex/",l.all[i],sep=""))
nm[i] <- strsplit(l.all[i],".",fixed=TRUE)[[1]][1]
image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i])
points(S)
points(S[which(is.na(X[,14])),],col="green")

x11()
image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i],xlim=c(-10^5,0),ylim=c(6400000,6700000))
points(S)
points(S[which(is.na(X[,14])),],col="green")

x11()
image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i],xlim=c(-4*10^4,-1*10^4),ylim=c(6570000,6600000))
points(S)
points(S[which(is.na(X[,14])),],col="green")



i=17

x11()
a[[i]] <- nc4.matrix(paste("./inputs/til_Alex/",l.all[i],sep=""))
nm[i] <- strsplit(l.all[i],".",fixed=TRUE)[[1]][1]
image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i])
points(S)
x11()
image(a[[i]]$x,a[[i]]$y,a[[i]]$z,main=nm[i],xlim=c(-10^5,0),ylim=c(6400000,6700000))
points(S)




data <- cbind(x=rnorm(6.9e2), y=rnorm(6.9e2))
labels <- rnorm(6.9e2)
query <- cbind(x=rnorm(5e1), y=rnorm(5e1))

library(FNN)
get.nn <- function(data, labels, query) {
  nns <- get.knnx(data, query, k=1)
  labels[nns$nn.index]
}
(get.nn(data, labels, query))

