rm(list=ls())
library(SpatGEVBMA)

data_wd="/nr/project/stat/ClimDesign/WP3/"

output.path <- "/nr/project/stat/ClimDesign/WP3/Res/Gumbel_latlon/"
covariates.folder <- "/nr/project/stat/ClimDesign/WP3/Data/toAlex/Cov/"
station.annualMax.file <- "/nr/project/stat/ClimDesign/WP3/Data/toAlex/Obs/AM_final.xlsx"
return.period <- c(2,5,10,20,25,50,100,200)
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

# 10 min - sheet 1:
#station.annualMax.sheet <- 1
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_onlyPlu.txt"
#output.folder.name <- "SpatGEV.res.10min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#		    fixed.xi = fixed.xi)

# 15 min - sheet 2:
#station.annualMax.sheet <- 2
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_onlyPlu.txt"
#output.folder.name <- "SpatGEV.res.15min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#  		    fixed.xi = fixed.xi)

# 20 min - sheet 3:
#station.annualMax.sheet <- 3
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_onlyPlu.txt"
#output.folder.name <- "SpatGEV.res.20min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                  station.annualMax.file = station.annualMax.file,
#                  station.annualMax.sheet = station.annualMax.sheet,
#                  station.locations.file = station.locations.file,
#                  output.path = output.path,
#                  output.folder.name = output.folder.name,
#                  return.period = return.period,
#                  post.quantiles = post.quantiles,
#                  transform.output = transform.output,  
#                  show.uncertainty = show.uncertainty,
#                  coordinate.type = coordinate.type,
#                  table.format = table.format,
#                  mcmc.reps = mcmc.reps,
#                  burn.in = burn.in,
#                  cores = cores,
#                  annualMax.name = annualMax.name,
#                  create.tempfiles = create.tempfiles,
#                  keep.temp.files = keep.temp.files,
#                  save.all.output = save.all.output,
#                  testing = testing,
#		    fixed.xi = fixed.xi)

# 30 min - sheet 4:
#station.annualMax.sheet <- 4
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_onlyPlu.txt"
#output.folder.name <- "SpatGEV.res.30min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#		    fixed.xi = fixed.xi)


# 45 min - sheet 5:
#station.annualMax.sheet <- 5
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_onlyPlu.txt"
#output.folder.name <- "SpatGEV.res.45min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#		   fixed.xi = fixed.xi)

# 60 min - sheet 6:
#station.annualMax.sheet <- 6
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_plu_and_geo.txt"
#output.folder.name <- "SpatGEV.res.60min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#fixed.xi = fixed.xi)

# 90 min - sheet 7:
#station.annualMax.sheet <- 7
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_onlyPlu.txt"
#output.folder.name <- "SpatGEV.res.90min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                  station.annualMax.file = station.annualMax.file,
#                  station.annualMax.sheet = station.annualMax.sheet,
#                  station.locations.file = station.locations.file,
#                  output.path = output.path,
#                  output.folder.name = output.folder.name,
#                  return.period = return.period,
#                  post.quantiles = post.quantiles,
#                  transform.output = transform.output,  
#                  show.uncertainty = show.uncertainty,
#                  coordinate.type = coordinate.type,
#                  table.format = table.format,
#                  mcmc.reps = mcmc.reps,
#                  burn.in = burn.in,
#                  cores = cores,
#                  annualMax.name = annualMax.name,
#                  create.tempfiles = create.tempfiles,
#                  keep.temp.files = keep.temp.files,
#                  save.all.output = save.all.output,
#                  testing = testing,
#fixed.xi = fixed.xi)

# 120 min - sheet 8:
#station.annualMax.sheet <- 8
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_plu_and_geo.txt"
#output.folder.name <- "SpatGEV.res.120min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                  station.annualMax.file = station.annualMax.file,
#                  station.annualMax.sheet = station.annualMax.sheet,
#                  station.locations.file = station.locations.file,
#                  output.path = output.path,
#                  output.folder.name = output.folder.name,
#                  return.period = return.period,
#                  post.quantiles = post.quantiles,
#                  transform.output = transform.output,  
#                  show.uncertainty = show.uncertainty,
#                  coordinate.type = coordinate.type,
#                  table.format = table.format,
#                  mcmc.reps = mcmc.reps,
#                  burn.in = burn.in,
#                  cores = cores,
#                 annualMax.name = annualMax.name,
#                 create.tempfiles = create.tempfiles,
#                 keep.temp.files = keep.temp.files,
#                 save.all.output = save.all.output,
#                 testing = testing,
#	   fixed.xi = fixed.xi)

# 180 min - sheet 9:
#station.annualMax.sheet <- 9
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_plu_and_geo.txt"
#output.folder.name <- "SpatGEV.res.180min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#		   fixed.xi = fixed.xi)

# 360 min - sheet 10:
#station.annualMax.sheet <- 10
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_plu_and_geo.txt"
#output.folder.name <- "SpatGEV.res.360min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#   fixed.xi = fixed.xi)

# 720 min - sheet 11:
#station.annualMax.sheet <- 11
#station.locations.file <- "~/lustreB/ExPrecFlood/NR/modelRunV2/Data/Meta/metadata_stations_plu_and_geo.txt"
#output.folder.name <- "SpatGEV.res.720min"

#SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
#                   station.annualMax.file = station.annualMax.file,
#                   station.annualMax.sheet = station.annualMax.sheet,
#                   station.locations.file = station.locations.file,
#                   output.path = output.path,
#                   output.folder.name = output.folder.name,
#                   return.period = return.period,
#                   post.quantiles = post.quantiles,
#                   transform.output = transform.output,  
#                   show.uncertainty = show.uncertainty,
#                   coordinate.type = coordinate.type,
#                   table.format = table.format,
#                   mcmc.reps = mcmc.reps,
#                   burn.in = burn.in,
#                   cores = cores,
#                   annualMax.name = annualMax.name,
#                   create.tempfiles = create.tempfiles,
#                   keep.temp.files = keep.temp.files,
#                   save.all.output = save.all.output,
#                   testing = testing,
#		   fixed.xi = fixed.xi)

# 1440 min - sheet 12:
station.annualMax.sheet <- 12
station.locations.file <- "/nr/project/stat/ClimDesign/WP3//Data/toAlex/Meta/metadata_stations_plu_and_geo.txt"
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
                   fixed.xi = fixed.xi)



