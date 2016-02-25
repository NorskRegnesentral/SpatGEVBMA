##  This script replicates the results from Dyrrdal et al.
rm(list = ls())

library(SpatialGEVBMA)

setwd("~/NR/SpatGEV/")

data(norway)
attach(norway)

Y <- Y.list
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


n.reps <- 2e5

R <- spatial.gev.bma(Y, X, S, n.reps, prior, print.every = 1e3)

tbl <- gev.process.results(R)

save(R, tbl, file="./output/gev.output.bma.RData")
