library(data.table)
setwd("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/CD_WP3/ExploreData")

#Params:
load("LocalGEVfits.RData")


metastation=data.table(read.table("/nr/project/stat/ClimDesign/WP3/Data/fromAnita/meta/metadata_new.txt",header=TRUE))
setnames(metastation,"Stnr","station")
metastation=metastation[!station=="63420...80"]
metastation=metastation[!station=="63420...132"]
metastation$station=as.numeric(metastation$station)


merge(metastation,DT,by="station")
