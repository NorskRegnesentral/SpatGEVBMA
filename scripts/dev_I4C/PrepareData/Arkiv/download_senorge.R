library(seNorge)
library(rgdal) 
library(data.table)
library(foreach)

setwd("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/")
AM60=fread(file="Data/AM60.csv")

preciploc=unique(AM60[,.(lon,lat,utmx,utmy,stid)])

to_get_coords=FALSE

if(to_get_coords){
  sN_coords=get_coord_key()
  
  xy <- SpatialPoints(cbind(Longitude=sN_coords$lon,Latitude=sN_coords$lat),proj4string = CRS("+proj=longlat")) 
  utm_sN=spTransform(xy, "+proj=utm +zone=33 +datum=WGS84") 
  utm_sN=utm_sN@coords
  
  radi=20000
  
  #-------Find seNorge locations to include--------------#
  ind_tokeep=c()
  statnr=c()
  for(j in 1:(dim(preciploc)[1])){
    stationloc=preciploc[j,.(utmx,utmy)]
    currkeep=which(((stationloc$utmx-utm_sN[,1])^2+(stationloc$utmy-utm_sN[,2])^2)<radi^2)
    ind_tokeep=c(ind_tokeep,currkeep)
    statnr=c(statnr,rep(preciploc$stid[j],length(currkeep)))
  }
  sN_coords_keep=data.table(sN_coords[ind_tokeep,])
  sN_unique=unique(sN_coords_keep)
  sN_coords_keep[,stid:=statnr]
  
  fwrite(sN_coords_keep,"Data/sN_covariates/sN_target_locations.csv")
  #---------------------------------------------------------#
}else{
  sN_coords_keep=fread("Data/sN_covariates/sN_target_locations.csv")
}
#Four seasons:
#juni,juli, august.
#sept,okt,nov.
#des,jan,feb.
#mars,apr,mai.


#----------------------------------
#https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
parallel::detectCores()
n.cores <- 8 #parallel::detectCores()-1

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)#PSOCK or FORK
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
#----------------------------------

startyear=1970;endyear=1974
tosave=list()
stationids=unique(sN_coords_keep$stid)

tosave=foreach(j =1:length(stationids))%dopar%{
  print(paste0("Location number: ",j))
  currcoords=sN_coords_keep[stid==stationids[j],]
  sNdata=load_sN(years=(startyear-1):endyear,locs=currcoords)
  sNmean=sNdata[,.("prec"=mean(prec),"st"=mean(st),"st_min"=mean(st_min),"st_max"=mean(st_max)),date]
  sNmean[,wetday:=(prec>=1)] #wetday if mean precip in area is over 1 mm.
  sNmean[,wetterday:=(prec>=5)] #wetday if mean precip in area is over 1 mm.
  sNmean[,verywetday:=(prec>=10)] #wetday if mean precip in area is over 1 mm.
  
  sNmean[,month:=month(date)]
  sNmean[,year:=year(date)]
  
  sN_summer=sNmean[month%in% c(6,7,8),.("precip_JJA"=sum(prec),"wetdays_JJA"=sum(wetday),"wetterdays_JJA"=sum(wetterday),"verywetday_JJA"=sum(verywetday),"temp_JJA"=mean(st,na.rm=TRUE),"max_temp_JJA"=mean(st_max),"min_temp_JJA"=mean(st_min)),.(year)]
  sN_fall=sNmean[month%in% c(9,10,11),.("precip_SON"=sum(prec),"wetdays_SON"=sum(wetday),"wetterdays_SON"=sum(wetterday),"verywetday_SON"=sum(verywetday),"temp_SON"=mean(st,na.rm=TRUE),"max_temp_SON"=mean(st_max),"min_temp_SON"=mean(st_min)),.(year)]
  sN_spring=sNmean[month%in% c(3,4,5),.("precip_MAM"=sum(prec),"wetdays_MAM"=sum(wetday),"wetterdays_MAM"=sum(wetterday),"verywetday_MAM"=sum(verywetday),"temp_MAM"=mean(st,na.rm=TRUE),"max_temp_MAM"=mean(st_max),"min_temp_MAM"=mean(st_min)),.(year)]
  
  sNmean[month==12,year:=year+1]
  sN_winter=sNmean[month%in% c(12,1,2),.("precip_DJF"=sum(prec),"wetdays_DJF"=sum(wetday),"wetterdays_DJF"=sum(wetterday),"verywetday_DJF"=sum(verywetday),"temp_DJF"=mean(st),"max_temp_DJF"=mean(st_max),"min_temp_DJF"=mean(st_min)),.(year)]

  sN_all=merge(sN_summer,sN_fall,by=c("year"),all.x=TRUE,all.y=TRUE)
  sN_all=merge(sN_all,sN_spring,by=c("year"),all.x=TRUE,all.y=TRUE)
  sN_all=merge(sN_all,sN_winter,by=c("year"),all.x=TRUE,all.y=TRUE)
  sN_all[,stid:=stationids[j]]
  
  sN_all=sN_all[year!=(startyear-1)]
  sN_all=sN_all[year!=(endyear+1)]
  
  sN_all=sN_all[is.na(min_temp_DJF)==FALSE,]
  
  
  return(sN_all)
}


tosave=rbindlist(tosave)
tosave=tosave[!(is.na(precip_DJF)==TRUE)]

fwrite(tosave,paste0("Data/sN_covariates/seasonal_weather_",min(tosave$year),"_",max(tosave$year),".csv"))

parallel::stopCluster(cl = my.cluster)
