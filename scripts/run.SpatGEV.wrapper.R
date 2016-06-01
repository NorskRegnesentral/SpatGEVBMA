rm(list=ls())

library(SpatialGEVBMA)

### Quick Test

covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.annualMax.sheet <- 3
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
output.folder.name <- "output_test5"
return.period <- c(20,100,200)
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 10^5 # Should at least be 10^5
burn.in <- 10^4
cores <- 10 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- FALSE
testing <- 500


SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.annualMax.file = station.annualMax.file,
                station.annualMax.sheet = station.annualMax.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
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
                testing = testing)

quit(save = "no")


### The full run 1


# Sheet 1

rm(list=ls())

library(SpatialGEVBMA)

station.annualMax.sheet <- 1 
return.period <- 20

covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 2*10^5 
burn.in <-4*10^4
cores <- 10 # 20 
annualMax.name <- NULL 
create.tempfiles <- TRUE
keep.temp.files <- TRUE
save.all.output <- TRUE
testing <- FALSE 

output.folder.name <- paste("output_dur_",station.annualMax.sheet,"_per_",return.period,sep="")

SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.annualMax.file = station.annualMax.file,
                station.annualMax.sheet = station.annualMax.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
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
                testing = testing)

## Sheet 2

rm(list=ls())

library(SpatialGEVBMA)


station.annualMax.sheet <- 2 
return.period <- 20
covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 2*10^5 
burn.in <-4*10^4
cores <- 10 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- TRUE
testing <- FALSE 

output.folder.name <- paste("output_dur_",station.annualMax.sheet,"_per_",return.period,sep="")

SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.annualMax.file = station.annualMax.file,
                station.annualMax.sheet = station.annualMax.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
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
                testing = testing)






## Sheet 4

rm(list=ls())

library(SpatialGEVBMA)


station.annualMax.sheet <- 4
return.period <- 20
covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 2*10^5 
burn.in <-4*10^4
cores <- 10 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- TRUE
testing <- FALSE 

output.folder.name <- paste("output_dur_",station.annualMax.sheet,"_per_",return.period,sep="")

SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.annualMax.file = station.annualMax.file,
                station.annualMax.sheet = station.annualMax.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
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
                testing = testing)


## Sheet 5

rm(list=ls())

library(SpatialGEVBMA)


station.annualMax.sheet <- 5
return.period <- 20
covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 2*10^5 
burn.in <-4*10^4
cores <- 10 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- TRUE
testing <- FALSE 

output.folder.name <- paste("output_dur_",station.annualMax.sheet,"_per_",return.period,sep="")

SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.annualMax.file = station.annualMax.file,
                station.annualMax.sheet = station.annualMax.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
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
                testing = testing)


## Sheet 6

rm(list=ls())

library(SpatialGEVBMA)


station.annualMax.sheet <- 6
return.period <- 20
covariates.folder <- "~/NR/SpatGEV/inputs/nc_files_used" 
station.annualMax.file <- "~/NR/SpatGEV/inputs/station_data/AM_allDurations.xlsx"
station.locations.file <- "~/NR/SpatGEV/inputs/station_data/metadata_stations_1hour.txt"
output.path <- "~/NR/SpatGEV"
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 2*10^5 
burn.in <-4*10^4
cores <- 10 # 20 
annualMax.name <- NULL 
create.tempfiles <- FALSE
keep.temp.files <- FALSE
save.all.output <- TRUE
testing <- FALSE 

output.folder.name <- paste("output_dur_",station.annualMax.sheet,"_per_",return.period,sep="")

SpatGEV.wrapper(covariates.folder = covariates.folder, 
                station.annualMax.file = station.annualMax.file,
                station.annualMax.sheet = station.annualMax.sheet, 
                station.locations.file = station.locations.file,
                output.path = output.path,
                output.folder.name = output.folder.name,
                return.period = return.period,
                post.quantiles = post.quantiles,
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
                testing = testing)











