library(SpatGEVBMA)
library(data.table)
library(devtools)

setwd("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/")
source("R/gev.R")
source("R/temporal_spatgev.R")
source("R/RcppExports.R")


amax_data=fread(file="scripts/dev_I4C/Data/AM1440_cov.csv")[,.(lon,lat,year,masl,stid,
                                                             wetterdays_MAM,wetterdays_JJA,wetterdays_SON,wetterdays_annual,
                                                             precip_MAM,precip_JJA,precip_SON,precip_annual,
                                                             temp_MAM,temp_JJA,temp_SON,temp_annual,y)]


amax_data=amax_data[lon<7&lat<62 & lat>58] #select locations in the Bergen area.

spatgev_data=make_temporal_spatgev_data(amax_data,TRUE)
mcmc_res=spatial.gev.bma(Y.list=spatgev_data$Y,X.all=spatgev_data$X,S=spatgev_data$S,n.reps=50000,
                         temporal=TRUE,print.every=100,nonspatial = FALSE)


save(mcmc_res,file="/nr/project/stat/Impetus4Change/Res/mcmc_temporal_spatial50000_daily.Rdata")

#Go to folder where the script is.
#export OPENBLAS_NUM_THREADS=1
#R CMD BATCH daily_test.R &
#cat daily_test.Rout
#cd "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/Test/"
#129320

