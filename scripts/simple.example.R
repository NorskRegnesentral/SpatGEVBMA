## A Simple example script that confirms that xi.constrain works.
rm(list = ls())

library(SpatGEVBMA)

data(norway)

mod2 = spatial.gev.bma(norway$Y.list, norway$X, norway$S, n.reps = 1e2, xi.constrain = c(0,.25))

zz = gev.impute(mod2, X.drop = norway$X[1:2,],S.drop = norway$S[1,] + rnorm(2,0,.1), xi.constrain = c(0,.25), return.param = TRUE)

print(all( zz$XI >= 0 & zz$XI <= 0.25))
