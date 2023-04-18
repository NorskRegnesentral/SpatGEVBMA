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
#source("R/gev.R")
#source("R/temporal_spatgev.R")
#-----------------------------------------------------#

#Upload hourly data (change 60 to 1440 to get daily data):
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")

#Select some covariates:
amax_data=amax_data[,.(year,masl,stid,lon,lat,
                       wetterdays_JJA_static_clim,wetterdays_annual_static_clim,
                       precip_JJA_static_clim,precip_annual_static_clim,
                       temp_JJA_static_clim,temp_annual_static_clim,y)]

#Select locations in the Bergen area:
amax_data=amax_data[lon<7&lat<62 & lat>58] 
amax_data=amax_data[year>=1972 & year<=2015 ]

uniquelocs=unique(amax_data[,stid])
Y=list()
for(j in 1:length(uniquelocs)){
  Y[[j]]=amax_data[stid==uniquelocs[j],y]
}

X=unique(amax_data[,.(masl,lon,lat,wetterdays_JJA_static_clim,wetterdays_annual_static_clim,
               precip_JJA_static_clim,precip_annual_static_clim,
               temp_JJA_static_clim,temp_annual_static_clim)])

X=cbind(1,scale(X))

S=unique(amax_data[,.(lon,lat)])
S=scale(S)


n.reps=100000 #mcmc iterations.
nonspatial=TRUE #Set nonspatial=TRUE to make the random effect iid.

#run spatgev with all years:
mcmc_all=spatial.gev.bma(Y.list=Y,X.all=X,S=S,n.reps=n.reps,
                         temporal=FALSE,print.every=1000,nonspatial = nonspatial)


save(mcmc_all,file=paste0("/nr/project/stat/Impetus4Change/Res/cv_spatial/nonspat",as.numeric(nonspatial),"_allyears.Rdata"))


#Then CV where we leave out each of these combinations:
cv_years=c(1980,1985,1990,1995,2000,2005,2010,2015)
for(j in 1:length(cv_years)){
  curr_amax=copy(amax_data)
  to_remove=cv_years[j]:(cv_years[j]+2)
  curr_amax=curr_amax[!(year%in%to_remove),]
  
  Y=list()
  for(j in 1:length(uniquelocs)){
    Y[[j]]=curr_amax[stid==uniquelocs[j],y]
  }
  
  mcmc_cv=spatial.gev.bma(Y.list=Y,X.all=X,S=S,n.reps=n.reps,
                          temporal=FALSE,print.every=1000,nonspatial = nonspatial)
  
  save(mcmc_cv,file=paste0("/nr/project/stat/Impetus4Change/Res/cv_spatial/nonspat",as.numeric(nonspatial),"_cv_",to_remove[1],"_",tail(to_remove,1),".Rdata"))
}





