library(SpatGEVBMA)
library(sp)

#----------------Prepare metadata-----------------------------------#
metadata_raw=fread("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/meta/metadata_raw_new.txt")
metadata_old=data.table(read.table("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/meta/metadata_stations_old.txt",header=TRUE))

metadata_all=merge(metadata_raw,metadata_old,by=c("Stnr"),all.x=TRUE)
metadata_all[is.na(Lon.y)==TRUE,Lon:=Lon.x]
metadata_all[is.na(Lon.y)==FALSE,Lon:=Lon.y]
metadata_all[is.na(Lat.y)==TRUE,Lat:=Lat.x]
metadata_all[is.na(Lat.y)==FALSE,Lat:=Lat.y]
metadata_all[is.na(Elev.y)==TRUE,Elev:=Elev.x]
metadata_all[is.na(Elev.y)==FALSE,Elev:=Elev.y]


metadata_all=metadata_all[,.(Stnr,X,Y,Lon,Lat,Elev)]

to_impute=metadata_all[is.na(X)==TRUE]
not_impute=metadata_all[is.na(X)==FALSE]


xy <- SpatialPoints(cbind(Longitude=to_impute$Lon,Latitude=to_impute$Lat),proj4string = CRS("+proj=longlat"))
res=spTransform(xy, "+proj=utm +zone=33 +datum=WGS84")

to_impute[,X:=round(res@coords[,1],0)]
to_impute[,Y:=round(res@coords[,1],0)]

all_data=rbind(to_impute,not_impute)
all_data=all_data[order(Stnr)]

all_data=all_data[,.(Stnr,Lon,Lat,X,Y,Elev)]
all_data=all_data[-c(28,63)]
write.table(all_data,"/nr/project/stat/ClimDesign/WP3/Data/fromOskar/meta/metadata_new.txt",row.names=FALSE)
#---------------------------------------------------------------------------------------------#


#----------------Prepare AM data-----------------------------------#
allsheets=list()
all_stations=unique(all_data$Stnr,2,100)


k=1
alldurs=c(10,15,20,30,45,60,90,120,180,360,720,1440)
for(dur in alldurs){
  curr_stations=all_stations
  current=fread(paste0("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/raw/obs/AM_",dur,"_min_all_plu_act_and_end_06_10_2020_v02.txt"),dec=",")
  
  datamat=matrix(NA,nrow=2020,ncol=length(all_stations))
  for(j in 1:length(all_stations)){
    obs=current[V1==all_stations[j],V4]
    years=current[V1==all_stations[j],V2]
    
    if(length(obs)!=0){
      datamat[years,j]=obs
    }else{
      print
    }
  }
  
  datamat=datamat[1970:2020,]
  datamat=cbind(1970:2020,datamat)
  
  colnames(datamat)=c("YEAR",all_stations)

  allsheets[[k]]=datamat
  k=k+1
}

names(allsheets)=paste0("AM_",alldurs,"min")

library("xlsx")
for(j in 1:length(allsheets)){
  write.xlsx(allsheets[[j]], file="/nr/project/stat/ClimDesign/WP3/Data/fromOskar/obs/AM_final.xlsx", 
             sheetName = names(allsheets)[j], 
             col.names = TRUE, row.names = FALSE, append = TRUE,showNA=FALSE)
}

#station.annualMax.file <- "/nr/project/stat/ClimDesign/WP3/Data/toAlex/Obs/AM_final.xlsx"
#fileYData <- XLConnect::loadWorkbook(station.annualMax.file)
#allYData <- suppressWarnings(XLConnect::readWorksheet(fileYData, sheet=station.annualMax.sheet))  # Ignore warnings
