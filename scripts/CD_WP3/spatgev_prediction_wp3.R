rm(list=ls())
library(ncdf4);library(Thermimage)
library(fields)
library(SpatGEVBMA)

#------------------------------------------------------------------------------------#
clim_years=2071:2100 #years we want to make predictions for. Possible values: 2021:2050, 2031:2060, 2041:2070, 2051:2080, 2061:2090, 2061:2100.
rcpnum=45 #rcp 26 or 45.
duration=1440 #duration.
#------------------------------------------------------------------------------------#

data_wd="/nr/project/stat/ClimDesign/WP3/"

# Take out parameter values from this file (mcmc run for rcpnum and duration):
load(paste0("/nr/project/stat/ClimDesign/WP3/Res/rcp",rcpnum,"/res_",duration,"min_rcp",rcpnum,"/mcmc.RData"))
mcmc.res=R0

#Take out covariates from this file (rcpnum and according to the first clim_year):
covariates.folder <- paste0("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/prediction/rcp",rcpnum,"/",clim_years[1],"/")


#Where to save the results:
output.path <- paste0("/nr/project/stat/ClimDesign/WP3/Res/rcp",rcpnum,"/prediction/",clim_years[1],"/")
output.folder.name <- paste0("res_",duration,"min","_rcp",rcpnum)



return.period = c(2,5,10,20,25,50,100,200)
post.quantiles = c(0.025,0.5,0.975)
xi.constrain = c(-Inf,Inf)
show.uncertainty=TRUE
create.tempfiles=TRUE
save.all.output=TRUE
annualMax.name <- duration
keep.temp.files=TRUE
fixed.xi=NULL
xi.constrain = c(-Inf,Inf)
testing=FALSE
burn.in <- 4*10^4
coordinate.type="XY"
transform.output = "UTM_33_to_LatLon"

SpatGEVBMA.wrapper.prediction(mcmc.res, #results file from .inference function.
                              covariates.folder, # Path to folder with covariate files in netcdf-format (see above) 
                              output.path =output.path,  # Path to the where the result folder should be stored
                              output.folder.name = output.folder.name,  # Name of result folder
                              return.period =return.period,  # Return period to impute results for (single number or a vector of numbers)
                              post.quantiles = post.quantiles,  # Vector of quantiles for which the posterior should be evaluated
                              show.uncertainty = show.uncertainty,  # Logical indicating whether an IQR uncertainty plot should also be provided
                              coordinate.type = coordinate.type, # Character indicating the type/name of coordinate system being used, either "XY" or "LatLon" (see above)
                              transform.output = transform.output, # Character specifying whether and how the output should be transformed. NULL corresponds to no transformation. "UTM_QQ_to_LatLon" transforms from UTM QQ (insert number) to LatLon
                              burn.in = burn.in, # The length of the initial burn-in period which is removed
                              cores = 16, # The number of cores on the computer used for the imputation. Using detectCores()-1 is good for running on a laptop.
                              annualMax.name = annualMax.name, # Name of annualMax data used in output plots and netcdf files. If NULL, then the name of the specified sheet is used.
                              create.tempfiles = create.tempfiles, # Logical indicating whether temporary files should be saved in a Temp folder to perform debugging and check intermediate variables/results if the function crashes
                              keep.temp.files = keep.temp.files, # Logical indicating whether the temporary files (if written) should be kept or deleted on function completion
                              save.all.output = save.all.output, # Logical indicating whether all R objects should be save to file upon function completion. Allocates approx 2.5 Gb for all of Norway.
                              testing = testing, # Variable indicating whether the run is a test or not. FALSE indicates no testing, a positive number indicates the number of locations being imputed
                              seed = 123, # The seed used in the mcmc computations
                              fixed.xi = fixed.xi,  # Where we want the shape parameter fixed
                              xi.constrain =xi.constrain)

#------------------------------------------------------------------------------------------------------------------#