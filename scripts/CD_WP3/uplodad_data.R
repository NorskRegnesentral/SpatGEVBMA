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
image.plot(mirror.matrix(tas[,,100]))
nc_close(nc)

#wetDays:
#nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.wetDays.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp45_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.wetDays.nc4")
wetDays=ncvar_get(nc,"wetDays")
image.plot(mirror.matrix(wetDays[,,90]))


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
nc=nc_open("https://thredds.met.no/thredds/dodsC/metusers/oskaral/ClimDesign/covariates/cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MAP.nc4")
nc=nc_open("cnrm-r1i1p1-aladin_rcp26_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MAP.nc4")

MSP=ncvar_get(nc,"MAP")
nc_close(nc)
