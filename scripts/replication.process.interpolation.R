rm(list = ls())

library(SpatGEVBMA)
library(ncdf4io)

setwd("~/NR/SpatGEV/")

load("./output/returns.RData")

load("./inputs/cov.RData")

Z.p <- matrix(unlist(l),ncol=3,byrow=TRUE)
w.na <- which(apply(is.na(cov),1,"any"))
Z <- matrix(NA,dim(cov)[1],5)
Z[-w.na,] <- cbind(cov[-w.na,2:1],Z.p)
Z[w.na,] <- cbind(cov[w.na,2:1], rbind(colMeans(Z.p),colMeans(Z.p)))

colnames(Z) <- c("lon","lat", .025,.5,.975)
rownames(Z) <- NULL


write.csv(Z, file="map.txt",row.names=FALSE)
 ##this all needs to be ncdf4.



