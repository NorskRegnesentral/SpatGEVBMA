library(ggplot2);library(data.table)
themeinfo=theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),legend.title = element_text(size=20),legend.text = element_text(size=20),plot.title = element_text(size = 20))


pca_data=function (variable, pressure_level = NA, descriptive_string = NULL,
                   out_dir = "/nr/project/stat/ClimateFutures/RenewableEnergy/SPMM/Volume/MastersData/")
{
  f_out = paste0(out_dir, variable, "_", ifelse(is.na(pressure_level),
                                                "", paste0(pressure_level, "hPa_")), ifelse(is.null(descriptive_string),
                                                                                            "", paste0(descriptive_string, "_")), "pca_analysis.RData")
  if (file.exists(f_out)) {
    load(f_out)
    l = list(dt_zonal = dt_zonal, dt_factor_loadings = dt_factor_loadings,
             dt_eigenvectors = dt_eigenvectors, X_mean = X_mean)
    return(l)
  }
  else {
    stop(paste0("File ", f_out, " not found"))
  }
}


#------------------------------------------------------------------------------------------------------------------------------#
amax_data60=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")[,.(lon,lat,year,masl,stid,
                                                               wetterdays_MAM,wetterdays_JJA,wetterdays_SON,wetterdays_annual,
                                                               precip_MAM,precip_JJA,precip_SON,precip_annual,
                                                               temp_MAM,temp_JJA,temp_SON,temp_annual,y)]

amax_data60=amax_data60[lon<7&lat<62 & lat>58] #select locations in the Bergen area.
ggplot(amax_data60,aes(x=year,y=y,group=as.factor(stid),col=as.factor(stid)))+geom_line()+themeinfo
#------------------------------------------------------------------------------------------------------------------------------#

amax_data1440=fread(file="scripts/dev_I4C/Data/AM1440_cov.csv")[,.(lon,lat,year,masl,stid,
                                                             wetterdays_MAM,wetterdays_JJA,wetterdays_SON,wetterdays_annual,
                                                             precip_MAM,precip_JJA,precip_SON,precip_annual,
                                                             temp_MAM,temp_JJA,temp_SON,temp_annual,y)]
amax_data1440=amax_data1440[lon<7&lat<62 & lat>58] #select locations in the Bergen area.


ggplot(amax_data1440,aes(x=year,y=y,group=as.factor(stid),col=as.factor(stid)))+geom_line()+themeinfo+xlim(c(1980,2020))
#------------------------------------------------------------------------------------------------------------------------------#

#50539, florida
amax_data1440[,.N,stid]
#47240, karmøy
amax_data[stid==50539]

#------------------------------------------------------------------------------------------------------------------------------#
load("/nr/project/stat/Impetus4Change/Data/Embeddings/specific_humidity_850hPa_atlantic_pca_analysis.RData")
dt_factor_loadings=dt_factor_loadings[,.(V1,V2,V3,V4,V5,V6,date,hour)]
dt_factor_loadings[,month:=month(date)]
dt_factor_loadings[,year:=year(date)]
dt_factor_loadings=dt_factor_loadings[month %in% c(1:12),]

dt_factor_loadings=dt_factor_loadings[,.("hum1"=mean(V1),"hum2"=mean(V2),"hum3"=mean(V3),"hum4"=mean(V4),"hum5"=mean(V5),"hum6"=mean(V6)),.(year)]


amax_data1440=merge(amax_data1440,dt_factor_loadings,"year")
amax_data60=merge(amax_data60,dt_factor_loadings,"year")


dat=pca_data("mean_sea_level_pressure_nao_standardized")
dt_factor_loadings_nao=dat$dt_factor_loadings

dt_factor_loadings_nao=dt_factor_loadings_nao[,.(V1,V2,V3,V4,V5,V6,date,hour)]
dt_factor_loadings_nao[,month:=month(date)]
dt_factor_loadings_nao[,year:=year(date)]
dt_factor_loadings_nao=dt_factor_loadings_nao[month %in% c(1:12),]

dt_factor_loadings_nao=dt_factor_loadings_nao[,.("nao1"=mean(V1),"nao2"=mean(V2),"nao3"=mean(V3),"nao4"=mean(V4),"nao5"=mean(V5),"nao6"=mean(V6)),.(year)]

amax_data1440=merge(amax_data1440,dt_factor_loadings_nao,"year")
amax_data60=merge(amax_data60,dt_factor_loadings_nao,"year")


#------------------------------------------------------------------------------------


#------------------------------------------#
#If we want to plot the whole country, not only bergen area:
if(0){
  amax_data60=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")[,.(lon,lat,year,masl,stid,
                                                                 wetterdays_MAM,wetterdays_JJA,wetterdays_SON,wetterdays_annual,
                                                                 precip_MAM,precip_JJA,precip_SON,precip_annual,
                                                                 temp_MAM,temp_JJA,temp_SON,temp_annual,y)]
  
  amax_data1440=fread(file="scripts/dev_I4C/Data/AM1440_cov.csv")[,.(lon,lat,year,masl,stid,
                                                                 wetterdays_MAM,wetterdays_JJA,wetterdays_SON,wetterdays_annual,
                                                                 precip_MAM,precip_JJA,precip_SON,precip_annual,
                                                                 temp_MAM,temp_JJA,temp_SON,temp_annual,y)]
  
  
  amax_data1440=merge(amax_data1440,dt_factor_loadings,"year")
  amax_data60=merge(amax_data60,dt_factor_loadings,"year")
  
  
  amax_data1440=merge(amax_data1440,dt_factor_loadings_nao,"year")
  amax_data60=merge(amax_data60,dt_factor_loadings_nao,"year")
}

ggplot(amax_data1440,aes(x=hum1,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data1440,aes(x=hum2,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data1440,aes(x=hum3,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data1440,aes(x=nao1,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data1440,aes(x=nao2,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data1440,aes(x=nao3,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)



ggplot(amax_data60,aes(x=hum1,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data60,aes(x=hum2,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data60,aes(x=hum3,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data60,aes(x=nao1,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data60,aes(x=nao2,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)

ggplot(amax_data60,aes(x=nao3,y=y))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x)




