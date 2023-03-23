library(data.table)
library(rgdal);library(sp)
setwd("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/")

#------------------------------------------------------------------------------------------#
#Annual max precip data with coordinates
AM60=fread("Data/AMnorway/AM_60_min_all_plu_act_and_end_06_10_2020_v02.txt",dec=",")
colnames(AM60)=c("stid","year","V3","y")
stationdata=fread("Data/AMnorway/Plu_stn_aktive_og avsluttet_06_10_2020_v02.csv")
setnames(stationdata,c("STNR","UTM_NORTH","UTM_EAST","NAME"),c("stid","utmy","utmx","name"))
stationdata=stationdata[,.(stid,utmx,utmy,name)]

AM60=merge(AM60,stationdata,"stid")

elev=fread("data/Data/AMnorway/elevlonlat.csv",dec=".")
setnames(elev,c("STNR","MASL","LON","LAT"),c("stid","masl","lon","lat"))
AM60=merge(AM60,unique(elev[,.(stid,masl,lon,lat)]),by=("stid"))
AM60[,n_obs:=.N,.(stid)]

#---How much data we have----#
n_stations=length(unique(AM60$stid))
n_obs=dim(AM60)[1]
minyear=min(AM60$year)
maxyear=max(AM60$year)

n_stations
minyear
maxyear

fwrite(AM60,file="data/AM60.csv")

#------------------------------------------------------------------------------------------#
