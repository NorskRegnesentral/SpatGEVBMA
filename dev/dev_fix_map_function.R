rm(list = ls())

library(SpatGEVBMA)

setwd("~/Desktop/toAlex/SpatGEV.res.1440min")

load("./Temp/print_maps_injection.RData")
j = 1

print_map(Q = Z.p[j,,],
          shortName = shortName,
          longName = longName,
          post.quantiles = post.quantiles,
          IQRLongName = IQRLongName,
          show.uncertainty = TRUE,
          ww.na = ww.na,
          n = n,
          output.x = output.x,
          output.y = output.y,
          output.name = output.name,
          filename.nc = filename.nc,
          filename.pdf = filename.pdf,
          main.quantile = main.quantile,
          main.iqr = main.iqr,
          nx = nx,
          ny = ny,
          dim.list = dim.list,
          all.post.quantiles = all.post.quantiles,
          transform.output = transform.output,
          original.image = original.image,
          XYGrid = XYGrid,
          coordinate.type = coordinate.type,
          S = S)

