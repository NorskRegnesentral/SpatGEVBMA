rm(list=ls())
library(ncdf4);library(Thermimage)
library(fields)

setwd("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/")

clim_years=1991:2020 #1967:2022, #1991:2020,#1971:2022 (clim 2)
type="inference"

#---------------------------------------------------------------------------------------#

#JJA temp:
nc=nc_open(paste0("raw/senorge/seNorge_1957-2022.JJAtemp.nc"))
lon=ncvar_get(nc,"lon")
lat=ncvar_get(nc,"lat")
Xc=ncvar_get(nc,"X")
Yc=ncvar_get(nc,"Y")
var=ncvar_get(nc,"JJAtemp")
years=1957:2022
nc_close(nc)

nc=nc_open(paste0("raw/rcp",45,"/cnrm-r1i1p1-aladin_rcp",45,"_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_1960-2100.JJAtemp.nc4"))
Xc=ncvar_get(nc,"Xc")
Yc=ncvar_get(nc,"Yc")
nc_close(nc)

toaverage=which(years%in%clim_years)
var_tokeep=var[,,toaverage]
var_mean=apply(var_tokeep,c(1,2),mean)

#-----Lage en ny ncdf--------------#
Xc <- ncdim_def( "X", "X", Xc)
Yc <- ncdim_def( "Y", "Y", Yc)

dimensions <- ncvar_def("tas", "Kelvin", dim=list(Xc,Yc), -1, 
                        longname="2m_temperature", prec="double")

ncnew <- nc_create( paste0(type,"/senorge/",clim_years[1],"/JJAtemp.nc"), dimensions)
ncvar_put( ncnew, dimensions, var_mean)
nc_close(ncnew)



#MAP
nc=nc_open(paste0("raw/senorge/seNorge_1957-2022.MAP.nc"))
lon=ncvar_get(nc,"lon")
lat=ncvar_get(nc,"lat")
var=ncvar_get(nc,"MAP")
years=1957:2022
nc_close(nc)

sek_in_year=60*60*24*365
var=var*sek_in_year

toaverage=which(years%in%clim_years)
var_tokeep=var[,,toaverage]
var_mean=apply(var_tokeep,c(1,2),mean)


#-----Lage en ny ncdf--------------#
Xc <- ncdim_def( "X", "X", Xc)
Yc <- ncdim_def( "Y", "Y", Yc)

dimensions <- ncvar_def("MAP", "mm", dim=list(Xc,Yc), -1, 
                        longname="Bias-Adjusted Precipitation", prec="double")

ncnew <- nc_create( paste0(type,"/senorge/",clim_years[1],"/MAP.nc"), dimensions)
ncvar_put( ncnew, dimensions, var_mean)
nc_close(ncnew)




#MSP
nc=nc_open(paste0("raw/senorge/seNorge_1957-2022.MSP.nc"))
lon=ncvar_get(nc,"lon")
lat=ncvar_get(nc,"lat")
var=ncvar_get(nc,"MSP")
years=1957:2022
nc_close(nc)

sek_in_summer=60*60*24*(30+31+30+31+31+30+31) #apr-oct.
var=var*sek_in_summer

toaverage=which(years%in%clim_years)
var_tokeep=var[,,toaverage]
var_mean=apply(var_tokeep,c(1,2),mean)

#-----Lage en ny ncdf--------------#
Xc <- ncdim_def( "X", "X", Xc)
Yc <- ncdim_def( "Y", "Y", Yc)

dimensions <- ncvar_def("MSP", "mm", dim=list(Xc,Yc), -1, 
                        longname="Bias-Adjusted Precipitation", prec="double")

ncnew <- nc_create( paste0(type,"/senorge/",clim_years[1],"/MSP.nc"), dimensions)
ncvar_put( ncnew, dimensions, var_mean)
nc_close(ncnew)




#MSP.JJA
nc=nc_open(paste0("raw/senorge/seNorge_1957-2022.MSPJJA.nc"))
lon=ncvar_get(nc,"lon")
lat=ncvar_get(nc,"lat")
var=ncvar_get(nc,"MSPJJA")
years=1957:2022
nc_close(nc)

sek_in_JJA=60*60*24*(30+31+31) #JJA
var=var*sek_in_JJA

toaverage=which(years%in%clim_years)
var_tokeep=var[,,toaverage]
var_mean=apply(var_tokeep,c(1,2),mean)

#-----Lage en ny ncdf--------------#
Xc <- ncdim_def( "X", "X", Xc)
Yc <- ncdim_def( "Y", "Y", Yc)

dimensions <- ncvar_def("MSPJJA", "mm", dim=list(Xc,Yc), -1, 
                        longname="Bias-Adjusted Precipitation", prec="double")

ncnew <- nc_create( paste0(type,"/senorge/",clim_years[1],"/MSPJJA.nc"), dimensions)
ncvar_put( ncnew, dimensions, var_mean)
nc_close(ncnew)




#wetDays
nc=nc_open(paste0("raw/senorge/seNorge_1957-2022.wetDays.nc"))
lon=ncvar_get(nc,"lon")
lat=ncvar_get(nc,"lat")
var=ncvar_get(nc,"wetDays")
years=1957:2022
nc_close(nc)

toaverage=which(years%in%clim_years)
var_tokeep=var[,,toaverage]
var_mean=apply(var_tokeep,c(1,2),mean)

#-----Lage en ny ncdf--------------#
Xc <- ncdim_def( "X", "X", Xc)
Yc <- ncdim_def( "Y", "Y", Yc)

dimensions <- ncvar_def("wetDays", "count", dim=list(Xc,Yc), -1, 
                        longname="number of days per year with at least 1 mm precipitatio", prec="double")

ncnew <- nc_create( paste0(type,"/senorge/",clim_years[1],"/wetDays.nc"), dimensions)
ncvar_put( ncnew, dimensions, var_mean)
nc_close(ncnew)


#---------------------------------------------------------------------------------------#
rcpnum=45
file.copy(from = paste0("raw/rcp",rcpnum,"/distSea.nc"),
          to = paste0(type,"/senorge/",clim_years[1],"/distSea.nc"))

file.copy(from = paste0("raw/rcp",rcpnum,"/elevation.nc"),
          to = paste0(type,"/senorge/",clim_years[1],"/elevation.nc"))

file.copy(from = paste0("raw/rcp",rcpnum,"/lon.nc"),
          to = paste0(type,"/senorge/",clim_years[1],"/lon.nc"))


file.copy(from = paste0("raw/rcp",rcpnum,"/lat.nc"),
          to = paste0(type,"/senorge/",clim_years[1],"/lat.nc"))
#---------------------------------------------------------------------------------------#

