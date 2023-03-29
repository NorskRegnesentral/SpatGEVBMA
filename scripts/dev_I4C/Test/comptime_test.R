rm(list = ls())

library(SpatGEVBMA)

#-----------------------------------------------------#
loc_thea = "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/"
if(file.exists(loc_thea)){
    setwd(loc_thea)
}else{
    setwd("~/pkg/SpatGEVBMA/")
}
source("R/gev.R")
source("R/temporal_spatgev.R")
#-----------------------------------------------------#

#Upload hourly data (change 60 to 1440 to get daily data):
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")

#Select some covariates:
amax_data=amax_data[,.(lon,lat,year,masl,stid,
                       wetterdays_MAM,wetterdays_JJA,wetterdays_SON,wetterdays_annual,
                       precip_MAM,precip_JJA,precip_SON,precip_annual,
                       temp_MAM,temp_JJA,temp_SON,temp_annual,y)]

#y is the response variable and should not be removed. Same with "year".


#Select locations in the Bergen area:
amax_data=amax_data[lon<7&lat<62 & lat>58] 
dim(amax_data)

#Data preparation:
spatgev_data=make_temporal_spatgev_data(amax_data,TRUE)

n.reps=50000 #mcmc iterations.
nonspatial=FALSE #Set nonspatial=TRUE to make the random effect iid.

#run spatgev:
mcmc_res=spatial.gev.bma(Y.list=spatgev_data$Y,X.all=spatgev_data$X,S=spatgev_data$S,n.reps=n.reps,
                         temporal=TRUE,print.every=100,nonspatial = nonspatial)


save(mcmc_res,file=paste0("/nr/project/stat/Impetus4Change/Res/test_n",n.reps,"_nonspat",as.numeric(nonspatial),".Rdata"))
