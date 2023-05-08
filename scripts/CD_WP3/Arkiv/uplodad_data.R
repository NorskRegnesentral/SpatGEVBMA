#https://thredds.met.no/thredds/catalog/metusers/oskaral/ClimDesign/covariates/catalog.html
library(ncdf4);library(Thermimage)
library(fields)

setwd("/nr/project/stat/ClimDesign/WP3/Data/")
#---rcp45----#
#JJA tmp:
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_1960-2100.JJAtemp.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_1960-2100.JJAtemp.nc4")
lon=ncvar_get(nc,"lon")
lat=ncvar_get(nc,"lon")
Xc=ncvar_get(nc,"Xc")
Yc=ncvar_get(nc,"Yc")
tas=ncvar_get(nc,"tas")
time=ncvar_get(nc,"time")
#projection_utm=ncvar_get(nc,"projection_utm")
image.plot(mirror.matrix(tas[,,1]))
nc_close(nc)

#wetDays:
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.wetDays.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.wetDays.nc4")
wetDays=ncvar_get(nc,"wetDays")
image.plot(mirror.matrix(wetDays[,,90]))
nc_close(nc)

#Bias adjusted precip (seasonal?):
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MSP.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MSP.nc4")
MSP=ncvar_get(nc,"MSP")
nc_close(nc)


#Bias adjusted precip (annual?):
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MAP.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MAP.nc4")
MSP=ncvar_get(nc,"MAP")
nc_close(nc)

#---rcp26---#
#JJA tmp:
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_1960-2100.JJAtemp.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_1960-2100.JJAtemp.nc4")

lon=ncvar_get(nc,"lon")
lat=ncvar_get(nc,"lon")
Xc=ncvar_get(nc,"Xc")
Yc=ncvar_get(nc,"Yc")
tas=ncvar_get(nc,"tas")
time=ncvar_get(nc,"time")
#projection_utm=ncvar_get(nc,"projection_utm")
image(mirror.matrix(tas[,,1]))
nc_close(nc)

#wetDays:
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.wetDays.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.wetDays.nc4")
wetDays=ncvar_get(nc,"wetDays")
nc_close(nc)


#Bias adjusted precip (seasonal?):
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MSP.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MSP.nc4")
MSP=ncvar_get(nc,"MSP")
nc_close(nc)


#Bias adjusted precip (annual?):
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MAP.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MAP.nc4")
MAP=ncvar_get(nc,"MAP")
nc_close(nc)

#-----------------------------------------------------------------------------#
#Look at Alex' data used before:
nc_alex=nc_open("toAlex/Cov/AMJJASOtemp.nc")
mean_temp=ncvar_get(nc_alex,"mean_temperature")
time=ncvar_get(nc_alex,"time_bnds")
X=ncvar_get(nc_alex,"X")
Y=ncvar_get(nc_alex,"Y")
nc_close(nc_alex)
dim(mean_temp)
#-----------------------------------------------------------------------------#


#-----Testfil fra Oskar--------#
nc_proj=nc_open("cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_1960-2100.JJAtemp.nc4")
tas=ncvar_get(nc_proj,"tas")
time=ncvar_get(nc_proj,"time")
Xc=ncvar_get(nc_proj,"Xc")
Yc=ncvar_get(nc_proj,"Yc")

#-----Lage en ny ncdf--------------#
Xc <- ncdim_def( "X", "X", Xc)
Yc <- ncdim_def( "Y", "Y", Yc)
time <- ncdim_def( "time", "time", time)

dimensions <- ncvar_def("tas", "Kelvin", dim=list(Xc,Yc,time), -1, 
                    longname="2m_temp", prec="double")

ncnew <- nc_create( "fromOskar/test.nc", dimensions)
ncvar_put( ncnew, dimensions, tas)
nc_close(ncnew)


#---Look at it-------#
look=nc_open("fromOskar/test.nc")
abc=ncvar_get(look,"tas")
xx=ncvar_get(look,"Xc")
yy=ncvar_get(look,"Yc")
tt=ncvar_get(look,"time")
image.plot(mirror.matrix(tas[,,2])-mirror.matrix(abc[,,2]))
#---------------------#

#Referanseperiode: 1991–2020
#Midten av århundret: 2041–2070
#Slutten av århundret: 2071–2100

#------------------------------------------------------#
#We want everything on this format:
nc_alex=nc_open("toAlex/Cov/AMJJASOtemp.nc")
tmp=ncvar_get(nc_alex,"mean_temperature")
XX=ncvar_get(nc_alex,"X")
YY=ncvar_get(nc_alex,"Y")
TT=ncvar_get(nc_alex,"time")
nc_close(nc_alex)
#------------------------------------------------------#

