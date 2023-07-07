library(SpatGEVBMA)
library(sp)

#----------------Prepare metadata-----------------------------------#
metadata_raw=fread("/nr/project/stat/ClimDesign/WP3/Data/fromAnita/raw/meta_stations.txt")


xy <- SpatialPoints(cbind(Longitude=metadata_raw$Lon,Latitude=metadata_raw$Lat),proj4string = CRS("+proj=longlat"))
res=spTransform(xy, "+proj=utm +zone=33 +datum=WGS84")

metadata_raw[,X:=round(res@coords[,1],0)]
metadata_raw[,Y:=round(res@coords[,2],0)]

all_data=metadata_raw[,.(Stnr,Lon,Lat,X,Y,Elev)]

write.table(all_data,"/nr/project/stat/ClimDesign/WP3/Data/fromAnita/meta/metadata_new.txt",row.names=FALSE)
#---------------------------------------------------------------------------------------------#
#just adding utmx and utmy.
