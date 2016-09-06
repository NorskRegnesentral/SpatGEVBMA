rm(list=ls())

library(SpatGEVBMA)

### Quick Test

covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.annualMax.sheet <- 3
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
output.folder.name <- "basic_100"
return.period <- c(100)
post.quantiles <- c(0.025,0.25,0.5,0.75,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
transform.output = NULL
table.format = "html"
mcmc.reps <- 300 # Should at least be 10^5
burn.in <- 100
cores <- 3 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- FALSE
testing <- 100
seed <- 123

SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
                   station.annualMax.file = station.annualMax.file,
                   station.annualMax.sheet = station.annualMax.sheet, 
                   station.locations.file = station.locations.file,
                   output.path = output.path,
                   output.folder.name = output.folder.name,
                   return.period = return.period,
                   post.quantiles = post.quantiles,
                   show.uncertainty = show.uncertainty,
                   coordinate.type = coordinate.type,
                   transform.output = transform.output,
                   table.format = table.format,
                   mcmc.reps = mcmc.reps,
                   burn.in = burn.in,
                   cores = cores,
                   annualMax.name = annualMax.name,
                   create.tempfiles = create.tempfiles,
                   keep.temp.files = keep.temp.files,
                   save.all.output = save.all.output,
                   testing = testing,
                   seed = seed)


rm(list=ls())

library(SpatGEVBMA)

### Quick Test

covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.annualMax.sheet <- 3
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
output.folder.name <- "basic_200"
return.period <- c(200)
post.quantiles <- c(0.025,0.25,0.5,0.75,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
transform.output = NULL
table.format = "html"
mcmc.reps <- 300 # Should at least be 10^5
burn.in <- 100
cores <- 3 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- FALSE
testing <- 100
seed <- 123

SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
                   station.annualMax.file = station.annualMax.file,
                   station.annualMax.sheet = station.annualMax.sheet, 
                   station.locations.file = station.locations.file,
                   output.path = output.path,
                   output.folder.name = output.folder.name,
                   return.period = return.period,
                   post.quantiles = post.quantiles,
                   show.uncertainty = show.uncertainty,
                   coordinate.type = coordinate.type,
                   transform.output = transform.output,
                   table.format = table.format,
                   mcmc.reps = mcmc.reps,
                   burn.in = burn.in,
                   cores = cores,
                   annualMax.name = annualMax.name,
                   create.tempfiles = create.tempfiles,
                   keep.temp.files = keep.temp.files,
                   save.all.output = save.all.output,
                   testing = testing,
                   seed = seed)



rm(list=ls())

library(SpatGEVBMA)

### Quick Test

covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.annualMax.sheet <- 3
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
output.folder.name <- "basic_100_200"
return.period <- c(100,200)
post.quantiles <- c(0.025,0.25,0.5,0.75,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
transform.output = NULL
table.format = "html"
mcmc.reps <- 300 # Should at least be 10^5
burn.in <- 100
cores <- 3 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- FALSE
testing <- 100
seed <- 123

SpatGEVBMA.wrapper(covariates.folder = covariates.folder,
                   station.annualMax.file = station.annualMax.file,
                   station.annualMax.sheet = station.annualMax.sheet, 
                   station.locations.file = station.locations.file,
                   output.path = output.path,
                   output.folder.name = output.folder.name,
                   return.period = return.period,
                   post.quantiles = post.quantiles,
                   show.uncertainty = show.uncertainty,
                   coordinate.type = coordinate.type,
                   transform.output = transform.output,
                   table.format = table.format,
                   mcmc.reps = mcmc.reps,
                   burn.in = burn.in,
                   cores = cores,
                   annualMax.name = annualMax.name,
                   create.tempfiles = create.tempfiles,
                   keep.temp.files = keep.temp.files,
                   save.all.output = save.all.output,
                   testing = testing,
                   seed = seed)


