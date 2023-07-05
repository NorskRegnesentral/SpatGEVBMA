library(data.table)
setwd("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/")

dur=1*60
sN_cov1=fread("Data/sN_covariates/seasonal_weather_1965_1969.csv")
sN_cov2=fread("Data/sN_covariates/seasonal_weather_1970_1979.csv")
sN_cov3=fread("Data/sN_covariates/seasonal_weather_1980_1989.csv")
sN_cov4=fread("Data/sN_covariates/seasonal_weather_1990_2021.csv")

sN_all=rbind(sN_cov1,sN_cov2,sN_cov3,sN_cov4)

fwrite(sN_all,"/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/Data/seasonal_sN_weather.csv")

#----------Make climatological variables---------#
sN_clim=list();k=1
startval=7;endval=7
for(j in 1967:2020){
  sN_curr=copy(sN_all[year%in%((j-startval):(j+endval))])
  sN_means=sN_curr[,.("precip_JJA_clim"=mean(precip_JJA),"temp_JJA_clim"=mean(temp_JJA),"wetdays_JJA_clim"=mean(wetdays_JJA),"wetterdays_JJA_clim"=mean(wetterdays_JJA),"verywetday_JJA_clim"=mean(verywetday_JJA),
             "wetdays_annual_clim"=mean(wetdays_annual),"wetterdays_annual_clim"=mean(wetterdays_annual),"verywetday_annual_clim"=mean(verywetday_annual),"precip_annual_clim"=mean(precip_annual),"temp_annual_clim"=mean(temp_annual)),
             .(stid)]
  sN_means[,year:=j,]
  
  sN_clim[[k]]=sN_means
  k=k+1
}
sN_clim=rbindlist(sN_clim)

sN_all=merge(sN_all,sN_clim,by=c("year","stid"))
sN_all[,climyears:=startval+endval+1]

#---------"Static climate, common for all"-----------------------------#
sN_all[,precip_JJA_static_clim:=mean(precip_JJA),.(stid)]
sN_all[,temp_JJA_static_clim:=mean(temp_JJA),.(stid)]
sN_all[,wetdays_JJA_static_clim:=mean(wetdays_JJA),.(stid)]
sN_all[,wetterdays_JJA_static_clim:=mean(wetterdays_JJA),.(stid)]
sN_all[,verywetday_JJA_static_clim:=mean(verywetday_JJA),.(stid)]

sN_all[,precip_annual_static_clim:=mean(precip_annual),.(stid)]
sN_all[,temp_annual_static_clim:=mean(temp_annual),.(stid)]
sN_all[,wetdays_annual_static_clim:=mean(wetdays_annual),.(stid)]
sN_all[,wetterdays_annual_static_clim:=mean(wetterdays_annual),.(stid)]
sN_all[,verywetday_annual_static_clim:=mean(verywetday_annual),.(stid)]

#------------------------------------------------------------------------#


amax_data=fread(file=paste0("Data/AM",dur,".csv"))[,.(lon,lat,year,masl,stid,y)]

amax_data=merge(amax_data,sN_all,by=c("stid","year"))


fwrite(amax_data,paste0("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/Data/AM",dur,"_cov.csv"))
