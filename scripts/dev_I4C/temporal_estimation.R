rm(list = ls())
library(SpatGEVBMA)
library(data.table)

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
amax_data=amax_data[,.(year,masl,stid,lon,lat,
                       wetterdays_JJA,wetterdays_annual,
                       precip_JJA,precip_annual,
                       temp_JJA,temp_annual,y)]

#Select locations in the Bergen area:
amax_data=amax_data[lon<7&lat<62 & lat>58] 

subset=TRUE
if(subset==TRUE){
  amax_data=amax_data[year>=1972 & year<=2015 ]
}

#Data preparation:
spatgev_data=make_temporal_spatgev_data(amax_data,TRUE)

n.reps=100000 #mcmc iterations.
nonspatial=TRUE #Set nonspatial=TRUE to make the random effect iid.

#run spatgev with all years:
mcmc_all=spatial.gev.bma(Y.list=spatgev_data$Y,X.all=spatgev_data$X,S=spatgev_data$S,n.reps=n.reps,
                         temporal=TRUE,print.every=1000,nonspatial = nonspatial)


if(subset==FALSE){
  save(mcmc_all,file=paste0("/nr/project/stat/Impetus4Change/Res/cv_temporal/nonspat",as.numeric(nonspatial),"_allyears.Rdata"))
}else{
  save(mcmc_all,file=paste0("/nr/project/stat/Impetus4Change/Res/cv_temporal_subset/nonspat",as.numeric(nonspatial),"_allyears.Rdata"))
}

#Then CV where we leave out each of these combinations:
cv_years=c(1980,1985,1990,1995,2000,2005,2010,2015)
for(j in 1:length(cv_years)){
  curr_amax=copy(amax_data)
  to_remove=cv_years[j]:(cv_years[j]+2)
  
  curr_amax=curr_amax[!(year%in%to_remove),]
  
  spatgev_data=make_temporal_spatgev_data(curr_amax,TRUE)
  
  mcmc_cv=spatial.gev.bma(Y.list=spatgev_data$Y,X.all=spatgev_data$X,S=spatgev_data$S,n.reps=n.reps,
                           temporal=TRUE,print.every=1000,nonspatial = nonspatial)
  
  if(subset==FALSE){
    save(mcmc_cv,file=paste0("/nr/project/stat/Impetus4Change/Res/cv_temporal/nonspat",as.numeric(nonspatial),"_cv_",to_remove[1],"_",tail(to_remove,1),".Rdata"))
  }else{
    save(mcmc_cv,file=paste0("/nr/project/stat/Impetus4Change/Res/cv_temporal_subset/nonspat",as.numeric(nonspatial),"_cv_",to_remove[1],"_",tail(to_remove,1),".Rdata"))
    
  }
}





