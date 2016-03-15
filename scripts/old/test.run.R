rm(list = ls())

library("SpatialGEVBMA")

data(norway)
attach(norway)
Rprof()
a <- spatial.gev.bma(Y.list, X, S, 2e1)
Rprof(NULL)
summaryRprof()


