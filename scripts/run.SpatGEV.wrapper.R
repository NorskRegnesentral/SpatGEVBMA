
rm(list=ls())

library(SpatialGEVBMA)

### Quick Test

covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.returns.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.returns.sheet <- 1 
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
output.folder.name <- "output_test2"
return.period <- 20
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
mcmc.reps <- 5*10^3 # Should at least be 10^4
burn.in <- 10^2
cores <- 10 # 20 
returns.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- FALSE
testing <- 10000 


SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.returns.file = station.returns.file,
                station.returns.sheet = station.returns.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
                show.uncertainty = show.uncertainty,
                coordinate.type = coordinate.type,
                mcmc.reps = mcmc.reps,
                burn.in = burn.in,
                cores = cores,
                returns.name = returns.name,
                create.tempfiles = create.tempfiles,
                keep.temp.files = keep.temp.files,
                save.all.output = save.all.output,
                testing = testing)



### The full run 1


# Sheet 1

rm(list=ls())

library(SpatialGEVBMA)

station.returns.sheet <- 1 
return.period <- 20

covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.returns.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
mcmc.reps <- 2*10^5 
burn.in <-4*10^4
cores <- 10 # 20 
returns.name <- NULL 
create.tempfiles <- TRUE
keep.temp.files <- TRUE
save.all.output <- TRUE
testing <- FALSE 

output.folder.name <- paste("output_dur_",station.returns.sheet,"_per_",return.period,sep="")

SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.returns.file = station.returns.file,
                station.returns.sheet = station.returns.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
                show.uncertainty = show.uncertainty,
                mcmc.reps = mcmc.reps,
                burn.in = burn.in,
                cores = cores,
                returns.name = returns.name,
                create.tempfiles = create.tempfiles,
                keep.temp.files = keep.temp.files,
                save.all.output = save.all.output,
                testing = testing)

## Sheet 2

rm(list=ls())

library(SpatialGEVBMA)


station.returns.sheet <- 2 
return.period <- 20
covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.returns.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
mcmc.reps <- 2*10^5 
burn.in <-4*10^4
cores <- 10 # 20 
returns.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- TRUE
testing <- FALSE 

output.folder.name <- paste("output_dur_",station.returns.sheet,"_per_",return.period,sep="")

SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.returns.file = station.returns.file,
                station.returns.sheet = station.returns.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
                show.uncertainty = show.uncertainty,
                mcmc.reps = mcmc.reps,
                burn.in = burn.in,
                cores = cores,
                returns.name = returns.name,
                create.tempfiles = create.tempfiles,
                keep.temp.files = keep.temp.files,
                save.all.output = save.all.output,
                testing = testing)




