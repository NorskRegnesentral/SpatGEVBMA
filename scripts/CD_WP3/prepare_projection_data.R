rm(list=ls())
library(ncdf4);library(Thermimage)
library(fields)

setwd("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/")
#source("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/CD_WP3/prepare_projection_data.R")


#1991:2020, 2021:2050, 2031:2060, 2041:2070, 2051:2080, 2061:2090, 2071:2100.
clim_years=2041:2070 #1991:2020,#1971:2022 (clim 2)
rcpnum=45
type="prediction"

climlist=list()
climlist[[1]]=1991:2020
climlist[[2]]=2021:2050
climlist[[3]]=2031:2060
climlist[[4]]=2041:2070
climlist[[5]]=2051:2080
climlist[[6]]=2061:2090
climlist[[7]]=2071:2100

for(pp in 1:length(climlist)){
clim_years=climlist[[pp]]
#print(clim_years)
for(type in c("prediction","inference")){

  
  if(type=="inference"){
    #if(pp==1){ clim_years=1991:2020}
    if(pp==1){clim_years=1967:2022}
  }
  
  if(pp>1 & type=="inference"){
    print("NEXT") 
    next
    }
  
  for(rcpnum in c(45)){
    print(type)
    print(clim_years)
  print(rcpnum)
    #JJA temp:
    
    nc=nc_open(paste0("raw/rcp",rcpnum,"/cnrm-r1i1p1-aladin_rcp",rcpnum,"_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_1960-2100.JJAtemp.nc4"))
    lon=ncvar_get(nc,"lon")
    lat=ncvar_get(nc,"lat")
    Xc=ncvar_get(nc,"Xc")
    Yc=ncvar_get(nc,"Yc")
    var=ncvar_get(nc,"tas")
    years=1960:2100
    nc_close(nc)
    
    
    toaverage=which(years%in%clim_years)
    var_tokeep=var[,,toaverage]
    var_mean=apply(var_tokeep,c(1,2),mean)
    
    #-----Lage en ny ncdf--------------#
    Xc <- ncdim_def( "X", "X", Xc)
    Yc <- ncdim_def( "Y", "Y", Yc)
    
    dimensions <- ncvar_def("tas", "Kelvin", dim=list(Xc,Yc), -1, 
                            longname="2m_temperature", prec="double")
    
    ncnew <- nc_create( paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/JJAtemp.nc"), dimensions)
    ncvar_put( ncnew, dimensions, var_mean)
    nc_close(ncnew)
    
    #---------------------------------------------------------------------------------------#
    
    #MAP temp:
    nc=nc_open(paste0("raw/rcp",rcpnum,"/cnrm-r1i1p1-aladin_rcp",rcpnum,"_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MAP.nc4"))
    lon=ncvar_get(nc,"lon")
    lat=ncvar_get(nc,"lat")
    Xc=ncvar_get(nc,"Xc")
    Yc=ncvar_get(nc,"Yc")
    var=ncvar_get(nc,"MAP")
    years=1960:2100
    nc_close(nc)
    
    sek_in_year=60*60*24*365
    var=var*sek_in_year
    
    toaverage=which(years%in%clim_years)
    var_tokeep=var[,,toaverage]
    var_mean=apply(var_tokeep,c(1,2),mean)
    
    #-----Lage en ny ncdf--------------#
    Xc <- ncdim_def( "X", "X", Xc)
    Yc <- ncdim_def( "Y", "Y", Yc)
    
    dimensions <- ncvar_def("MAP", "mm/yr", dim=list(Xc,Yc), -1, 
                            longname="Bias-Adjusted Precipitation", prec="double")
    
    ncnew <- nc_create( paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/MAP.nc"), dimensions)
    ncvar_put( ncnew, dimensions, var_mean)
    nc_close(ncnew)
    
    
    #---------------------------------------------------------------------------------------#
    
    #MSP temp:
    nc=nc_open(paste0("raw/rcp",rcpnum,"/cnrm-r1i1p1-aladin_rcp",rcpnum,"_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MSP.nc4"))
    lon=ncvar_get(nc,"lon")
    lat=ncvar_get(nc,"lat")
    Xc=ncvar_get(nc,"Xc")
    Yc=ncvar_get(nc,"Yc")
    var=ncvar_get(nc,"MSP")
    years=1960:2100
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
    
    ncnew <- nc_create( paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/MSP.nc"), dimensions)
    ncvar_put( ncnew, dimensions, var_mean)
    nc_close(ncnew)
    
    
    
    #JJA precip:
    nc=nc_open(paste0("raw/rcp",rcpnum,"/cnrm-r1i1p1-aladin_rcp",rcpnum,"_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.MSPJJA.nc4"))
    lon=ncvar_get(nc,"lon")
    lat=ncvar_get(nc,"lat")
    Xc=ncvar_get(nc,"Xc")
    Yc=ncvar_get(nc,"Yc")
    var=ncvar_get(nc,"MSPJJA")
    years=1960:2100
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
    
    ncnew <- nc_create( paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/MSPJJA.nc"), dimensions)
    ncvar_put( ncnew, dimensions, var_mean)
    nc_close(ncnew)
    
    
    #---------------------------------------------------------------------------------------#
    #wetDays:
    nc=nc_open(paste0("raw/rcp",rcpnum,"/cnrm-r1i1p1-aladin_rcp",rcpnum,"_eqm-sn2018v2005_rawbc_norway_1km_pr_daily_1960-2100.wetDays.nc4"))
    lon=ncvar_get(nc,"lon")
    lat=ncvar_get(nc,"lat")
    Xc=ncvar_get(nc,"Xc")
    Yc=ncvar_get(nc,"Yc")
    var=ncvar_get(nc,"wetDays")
    years=1960:2100
    nc_close(nc)
    
    
    toaverage=which(years%in%clim_years)
    var_tokeep=var[,,toaverage]
    var_mean=apply(var_tokeep,c(1,2),mean)
    
    #-----Lage en ny ncdf--------------#
    Xc <- ncdim_def( "X", "X", Xc)
    Yc <- ncdim_def( "Y", "Y", Yc)
    
    dimensions <- ncvar_def("wetDays", "count", dim=list(Xc,Yc), -1, 
                            longname="number of days per year with at least 1 mm precipitation", prec="double")
    
    ncnew <- nc_create( paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/wetDays.nc"), dimensions)
    ncvar_put( ncnew, dimensions, var_mean)
    nc_close(ncnew)
    #---------------------------------------------------------------------------------------#
    
    #wetDays:
    file.copy(from = paste0("raw/rcp",rcpnum,"/distSea.nc"),
              to = paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/distSea.nc"))
    
    #---------------------------------------------------------------------------------------#
    file.copy(from = paste0("raw/rcp",rcpnum,"/elevation.nc"),
              to = paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/elevation.nc"))
    
    file.copy(from = paste0("raw/rcp",rcpnum,"/lon.nc"),
              to = paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/lon.nc"))
    
    
    file.copy(from = paste0("raw/rcp",rcpnum,"/lat.nc"),
              to = paste0(type,"/rcp",rcpnum,"/",clim_years[1],"/lat.nc"))
  
  #elev, distsea, lon and lat are taken from Alex' data (toAlex). The rest is new data from Oskar.
  }
  }
}