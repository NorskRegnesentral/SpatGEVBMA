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

dat=pca_data("mean_sea_level_pressure_nao_standardized")

length(unique(dat$dt_zonal$lat))
length(unique(dat$dt_zonal$date))

#factor loadings kan brukes som kovariat:
plot(dat$dt_zonal$zonal_mean[1:1000],type="l")
plot(dat$dt_factor_loadings$V1[1:1000], type="l")
lines(dat$dt_factor_loadings$V2[1:1000],col="red")
lines(dat$dt_factor_loadings$V3[1:1000],col="blue")
lines(dat$dt_factor_loadings$V4[1:1000],col="gray")
lines(dat$dt_factor_loadings$V5[1:1000],col="gray")


dim(dat$dt_eigenvectors)



plot(dat$dt_factor_loadings$V1[1:5000],type="l")

dat$dt_zonal
dat$X_mean