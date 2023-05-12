rm(list=ls())
library(ncdf4);library(Thermimage)
library(fields)
library(SpatGEVBMA)

#------------------------------------------------------------------------------------#

clim_years=1991:2020 #1991:2020,#1971:2022 (clim 2)
rcpnum=45
duration=10

#------------------------------------------------------------------------------------#

durtable=data.table(dur=c(10,15,20,30,45,60,90,120,180,360,720,1440),sheet=1:12)
station.annualMax.sheet <- durtable[dur==duration,sheet]


data_wd="/nr/project/stat/ClimDesign/WP3/"
output.folder.name <- paste0("res_",duration,"min","_rcp",rcpnum)

output.path <- paste0("/nr/project/stat/ClimDesign/WP3/Res/rcp",rcpnum,"/")
covariates.folder <- paste0("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/inference/rcp",rcpnum,"/",clim_years[1],"/")
station.annualMax.file <- "/nr/project/stat/ClimDesign/WP3/Data/fromOskar/obs/AM_final.xlsx"
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


station.locations.file <- "/nr/project/stat/ClimDesign/WP3//Data/fromOskar/meta/metadata_new.txt"

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

