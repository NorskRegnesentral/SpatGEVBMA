rm(list=ls())

library(SpatGEVBMA)
setwd("~/Desktop/toAlex/")
      
output.path <- "./Results/"
covariates.folder <- "./Data/Cov/"
station.annualMax.file <- "./Data/Obs/AM_final.xlsx"
return.period <- c(2,5,10,20,25,30,40,50,100,200)
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 2*10^5 
burn.in <- 4*10^4
cores <- 16 # 20
annualMax.name <- NULL
create.tempfiles <- TRUE
keep.temp.files <- FALSE
save.all.output <- TRUE
testing <- FALSE
transform.output = "UTM_33_to_LatLon"
fixed.xi = NULL
xi.constrain = c(-0.25,0.3)


# 1440 min - sheet 12:
station.annualMax.sheet <- 12
station.locations.file <- "./Data/Meta/metadata_stations_plu_and_geo.txt"
output.folder.name <- "SpatGEV.res.1440min"

SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
                   station.annualMax.file = station.annualMax.file,
                   station.annualMax.sheet = station.annualMax.sheet,
                   station.locations.file = station.locations.file,
                   output.path = output.path,
                   output.folder.name = output.folder.name,
                   return.period = return.period,
                   post.quantiles = post.quantiles,
                   transform.output = transform.output,  
                   show.uncertainty = show.uncertainty,
                   coordinate.type = coordinate.type,
                   table.format = table.format,
                   mcmc.reps = mcmc.reps,
                   burn.in = burn.in,
                   cores = cores,
                  annualMax.name = annualMax.name,
                  create.tempfiles = create.tempfiles,
                  keep.temp.files = keep.temp.files,
                  save.all.output = save.all.output,
                  testing = testing,
                  fixed.xi = fixed.xi, xi.constrain = xi.constrain)
