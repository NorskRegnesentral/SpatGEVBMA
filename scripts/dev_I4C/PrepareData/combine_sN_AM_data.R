library(data.table)
setwd("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/")

sN_cov1=fread("Data/sN_covariates/seasonal_weather_1965_1969.csv")
sN_cov2=fread("Data/sN_covariates/seasonal_weather_1970_1979.csv")
sN_cov3=fread("Data/sN_covariates/seasonal_weather_1980_1989.csv")
sN_cov4=fread("Data/sN_covariates/seasonal_weather_1990_2021.csv")

sN_all=rbind(sN_cov1,sN_cov2,sN_cov3,sN_cov4)

fwrite(sN_all,"/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/Data/seasonal_sN_weather.csv")

amax_data=fread(file="Data/AM60.csv")[,.(lon,lat,year,masl,stid,y)]

amax_data=merge(amax_data,sN_all,by=c("stid","year"))


fwrite(amax_data,"/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/Data/AM60_cov.csv")
