library(akima);library(lubridate);library(ggplot2)
findPCA=function(X,dt_coord,scale=FALSE){
  #--------------------------------------------------#
  if(0){#Example
    load(file="/nr/samba/user/roksvag/Prosjekter2021/Onions/data/weatherdat_ERA5.Rdata")
    n_eig=5
    dt_index = unique(weatherPCA[,.(year)])
    dt_coord = unique(weatherPCA[, .(lon, lat)])
    
    X1 = matrix(weatherPCA[, prec],
                nrow = dt_index[, .N],
                byrow = TRUE)
    
    res=findPCA(X1,dt_coord,scale=TRUE)
    U=res
    dt_U=res$dt_U
    coef=res$pca_coef
    colnames(coef)=c(paste0("PC",1:n_eig))
    rownames(coef)=dt_index$year
  }
  #--------------------------------------------------#
  
  mu = colMeans(X)
  if(scale==FALSE){
    sig=apply(X,2,sd)
  }else{
    sig=mu*0+1
  }
  
  X_center = X
  for(j in 1:dim(X)[2]){
    X_center[, j] = (X[, j] - mu[j])/sig[j]
  }
  
  Sigma = t(X_center) %*% X_center
  l_svd = svd(Sigma)
  U = l_svd$u[, 1:n_eig]
  
  #PCA coefficient by year:
  coef=X_center%*%U
  
  dim(X_center) 
  dim(U)
  dt_U = cbind(dt_coord, U, mu)
  return(list(U=U,dt_U=dt_U,pca_coef=coef))
}



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


#-----------------------------------------------#
load("/nr/project/stat/Impetus4Change/Data/Embeddings/specific_humidity_850hPa_atlantic_pca_analysis.RData")
coords=unique(dt_eigenvectors[,.(lon,lat)])
dim(dt_eigenvectors)


s1=interp(x=dt_eigenvectors$lon,y=dt_eigenvectors$lat,z=dt_eigenvectors$V1)
s2=interp(x=dt_eigenvectors$lon,y=dt_eigenvectors$lat,z=dt_eigenvectors$V2)
s3=interp(x=dt_eigenvectors$lon,y=dt_eigenvectors$lat,z=dt_eigenvectors$V3)
s4=interp(x=dt_eigenvectors$lon,y=dt_eigenvectors$lat,z=dt_eigenvectors$V4)
s5=interp(x=dt_eigenvectors$lon,y=dt_eigenvectors$lat,z=dt_eigenvectors$V5)
s6=interp(x=dt_eigenvectors$lon,y=dt_eigenvectors$lat,z=dt_eigenvectors$V6)

par(mfrow=c(2,3),mar=c(5,5,5,5))
image.plot(s1);title("V1");image.plot(s2);title("V2");image.plot(s3);title("V3");
image.plot(s4);title("V4");image.plot(s5);title("V5");image.plot(s6);title("V6");


#------------------------------------------------------------------#
load("/nr/project/stat/Impetus4Change/Data/Embeddings/specific_humidity_850hPa_atlantic_pca_analysis.RData")
dt_factor_loadings=dt_factor_loadings[,.(V1,V2,V3,V4,V5,V6,date,hour)]
dt_factor_loadings[,month:=month(date)]
dt_factor_loadings[,year:=year(date)]

dt_factor_loadings=dt_factor_loadings[,.(mean(V1),mean(V2),mean(V3),mean(V4),mean(V5),mean(V6)),.(month,year)]


dt_factor_loadings

ggplot(dt_factor_loadings[year %in% seq(1980,2020,by=4)],aes(x=month,y=V1,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings[year %in% seq(1980,2020,by=4)],aes(x=month,y=V2,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings[year %in% seq(1980,2020,by=4)],aes(x=month,y=V3,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings[year %in% seq(1980,2020,by=4)],aes(x=month,y=V4,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings[year %in% seq(1980,2020,by=4)],aes(x=month,y=V5,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings[year %in% seq(1980,2020,by=4)],aes(x=month,y=V6,group=as.factor(year),color=as.factor(year)))+geom_line()


#------------------------------------------------------------------#
load("/nr/project/stat/Impetus4Change/Data/Embeddings/specific_humidity_850hPa_atlantic_pca_analysis.RData")
dt_factor_loadings=dt_factor_loadings[,.(V1,V2,V3,V4,V5,V6,date,hour)]
dt_factor_loadings[,month:=month(date)]
dt_factor_loadings[,year:=year(date)]
dt_factor_loadings=dt_factor_loadings[month %in% c(7,8,9),]

dt_factor_loadings=dt_factor_loadings[,.(mean(V1),mean(V2),mean(V3),mean(V4),mean(V5),mean(V6)),.(year)]


ggplot(dt_factor_loadings,aes(x=year,y=V1))+geom_line()+geom_point(pch=21,cex=3,bg="orange")
ggplot(dt_factor_loadings,aes(x=year,y=V2))+geom_line()+geom_point(pch=21,cex=3,bg="orange")
ggplot(dt_factor_loadings,aes(x=year,y=V3))+geom_line()+geom_point(pch=21,cex=3,bg="orange")
ggplot(dt_factor_loadings,aes(x=year,y=V4))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")
ggplot(dt_factor_loadings,aes(x=year,y=V5))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")
ggplot(dt_factor_loadings,aes(x=year,y=V6))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")





#------------------------------------------------------------------#
dat=pca_data("mean_sea_level_pressure_nao_standardized")
dt_eigenvectors_nao=dat$dt_eigenvectors

coords=unique(dt_eigenvectors_nao[,.(lon,lat)])
dim(dt_eigenvectors_nao)


s1=interp(x=dt_eigenvectors_nao$lon,y=dt_eigenvectors_nao$lat,z=dt_eigenvectors_nao$V1)
s2=interp(x=dt_eigenvectors_nao$lon,y=dt_eigenvectors_nao$lat,z=dt_eigenvectors_nao$V2)
s3=interp(x=dt_eigenvectors_nao$lon,y=dt_eigenvectors_nao$lat,z=dt_eigenvectors_nao$V3)
s4=interp(x=dt_eigenvectors_nao$lon,y=dt_eigenvectors_nao$lat,z=dt_eigenvectors_nao$V4)
s5=interp(x=dt_eigenvectors_nao$lon,y=dt_eigenvectors_nao$lat,z=dt_eigenvectors_nao$V5)
s6=interp(x=dt_eigenvectors_nao$lon,y=dt_eigenvectors_nao$lat,z=dt_eigenvectors_nao$V6)

par(mfrow=c(2,3),mar=c(5,5,5,5))
image.plot(s1);title("V1");image.plot(s2);title("V2");image.plot(s3);title("V3");
image.plot(s4);title("V4");image.plot(s5);title("V5");image.plot(s6);title("V6");

dt_factor_loadings_nao=dat$dt_factor_loadings

dt_factor_loadings_nao=dt_factor_loadings_nao[,.(V1,V2,V3,V4,V5,V6,date,hour)]
dt_factor_loadings_nao=dt_factor_loadings_nao[,.(mean(V1),mean(V2),mean(V3),mean(V4),mean(V5),mean(V6)),.(date)]

dt_factor_loadings_nao[,month:=month(date)]
dt_factor_loadings_nao[,year:=year(date)]

dt_factor_loadings_nao=dt_factor_loadings_nao[,.(mean(V1),mean(V2),mean(V3),mean(V4),mean(V5),mean(V6)),.(month,year)]

dt_factor_loadings_nao

ggplot(dt_factor_loadings_nao[year %in% seq(1980,2020,by=4)],aes(x=month,y=V1,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings_nao[year %in% seq(1980,2020,by=4)],aes(x=month,y=V2,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings_nao[year %in% seq(1980,2020,by=4)],aes(x=month,y=V3,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings_nao[year %in% seq(1980,2020,by=4)],aes(x=month,y=V4,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings_nao[year %in% seq(1980,2020,by=4)],aes(x=month,y=V5,group=as.factor(year),color=as.factor(year)))+geom_line()
ggplot(dt_factor_loadings_nao[year %in% seq(1980,2020,by=4)],aes(x=month,y=V6,group=as.factor(year),color=as.factor(year)))+geom_line()



#-------------------
dat=pca_data("mean_sea_level_pressure_nao_standardized")
dt_factor_loadings_nao=dat$dt_factor_loadings

dt_factor_loadings_nao=dt_factor_loadings_nao[,.(V1,V2,V3,V4,V5,V6,date,hour)]
dt_factor_loadings_nao[,month:=month(date)]
dt_factor_loadings_nao[,year:=year(date)]
dt_factor_loadings_nao=dt_factor_loadings_nao[month %in% c(9,10,11),]

dt_factor_loadings_nao=dt_factor_loadings_nao[,.(mean(V1),mean(V2),mean(V3),mean(V4),mean(V5),mean(V6)),.(year)]


ggplot(dt_factor_loadings_nao,aes(x=year,y=V1))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")
ggplot(dt_factor_loadings_nao,aes(x=year,y=V2))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")
ggplot(dt_factor_loadings_nao,aes(x=year,y=V3))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")
ggplot(dt_factor_loadings_nao,aes(x=year,y=V4))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")
ggplot(dt_factor_loadings_nao,aes(x=year,y=V5))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")
ggplot(dt_factor_loadings_nao,aes(x=year,y=V6))+geom_line()+geom_point(pch=21,cex=3,bg="skyblue")

