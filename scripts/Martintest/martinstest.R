### This is Martins test file, essentially a copy of test.run
# Doing this from Rstudio -- sorry Alex :)

rm(list = ls())

library("SpatialGEVBMA")

data(norway)
attach(norway)
Rprof()
a <- spatial.gev.bma(Y.list, X, S, 2e2)
Rprof(NULL)
head(summaryRprof()$by.total,15)

# Hierarchy:
#spatial.gev.bma
  # gev.update
    # gev.update.hyper
      # gev.update.lambda
        # l.prime
        # l.double.prime
    # gev.update.theta
    # gev.update.tau.xi 
      # j.prime
      # j.double.prime



