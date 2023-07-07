library(INLA)
library(evgam)
library(evd)
library(data.table)
loc_thea = "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/"
setwd(loc_thea)
library(SpatGEVBMA)
source("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/R/bgev.R")
#source("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/bGev_tutorial/cv_bgev.R")

#Same as the other, but now we let the tail and spread depend on precip_JJA and temp_JJA.
#------------------------------------------------------------------------------------------------------------------#
p.alpha = 0.5; p.beta = 0.25
spread = 0.3
tail = 0.1

#Setting priors:
hyper.spread = list(initial = 1, fixed=FALSE, prior="loggamma",param=c(3,3))

tail.interval= c(0,0.5)
tail.intern= map.tail(tail,tail.interval,inverse=TRUE)

hyper.tail = list(initial = tail.intern,
                  prior = "pc.gevtail",
                  param = c(7, tail.interval),
                  fixed= FALSE)

hyper.tail = list(initial = if (tail == 0.0) -Inf else tail.intern,
                  prior = "pc.gevtail",
                  param = c(7, tail.interval),
                  fixed= if (tail == 0.0) TRUE else FALSE)

hyper.bgev = list(spread = hyper.spread,
                  tail = hyper.tail)


control.bgev = list(q.location = p.alpha,
                    q.spread = p.beta,
                    # quantile levels for the mixing part
                    q.mix= c(0.05, 0.30),
                    # the Beta(s1, s2) mixing distribution parameters.
                    # Hard-coded for the moment: s1=s2=5
                    beta.ab = 5)
#------------------------------------------------------------------------------------------------------------------#

#Upload hourly data (change 60 to 1440 to get daily data):
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")[,.(year,masl,stid,lon,lat,
                                                             wetterdays_JJA,
                                                             precip_JJA,
                                                             temp_JJA,wetterdays_JJA_static_clim,precip_JJA_static_clim, temp_JJA_static_clim,year,y)]



reslist=list()
allyears=sort(unique(amax_data$year))
p=1
for(thisyear in 1980:2020){
  print(thisyear)
  amax2=copy(amax_data[,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,wetterdays_JJA_static_clim,precip_JJA_static_clim, temp_JJA_static_clim,masl,year)])
  cvdata=amax2[year==thisyear]
  amax2=amax2[year!=thisyear,]
  
  amax2= data.table(scale(amax2[,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl,wetterdays_JJA_static_clim,precip_JJA_static_clim, temp_JJA_static_clim)]))
  amax2[,intercept:=1]
  
  y_meansd=amax_data[year!=thisyear,.(mean_y=mean(y),sd_y=sd(y)),]
  covmeans=colMeans(amax_data[year!=thisyear,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl,year,wetterdays_JJA_static_clim,precip_JJA_static_clim, temp_JJA_static_clim)])
  covsd=apply(amax_data[year!=thisyear,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl,year,wetterdays_JJA_static_clim,precip_JJA_static_clim, temp_JJA_static_clim)],2,sd)
  
  cvdata=cvdata[,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl,year,wetterdays_JJA_static_clim,precip_JJA_static_clim,temp_JJA_static_clim)]
  amax2=amax2[,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl,year,wetterdays_JJA_static_clim,precip_JJA_static_clim,temp_JJA_static_clim,intercept)]
  cvdata0=cvdata
  cvdata=(cvdata-matrix(rep(covmeans,each=dim(cvdata)[1]),byrow=FALSE,ncol=dim(cvdata)[2]))/matrix(rep(covsd,each=dim(cvdata)[1]),byrow=FALSE,ncol=dim(cvdata)[2])
  cvdata[,intercept:=1]
  cvdata[,y:=NA]
  cvdata[,year:=NULL]
  amax2[,year:=NULL]
  
  
  amax2=rbind(amax2,cvdata)
  
  lon=amax2$lon;lat=amax2$lat
  formula1 = inla.mdata(y,cbind(lon,lat,temp_JJA) ,cbind(lon,lat,temp_JJA)) ~ -1 + intercept + wetterdays_JJA+precip_JJA+temp_JJA+masl+lon+lat
  formula2 = inla.mdata(y,cbind(lon,lat,temp_JJA_static_clim),cbind(lon,lat,temp_JJA_static_clim)) ~ -1 + intercept + wetterdays_JJA_static_clim+precip_JJA_static_clim+temp_JJA_static_clim+masl+lon+lat
  
  inlares1 = inla(formula1,
                  family = "bgev",
                  data = amax2,
                  control.family = list(hyper = hyper.bgev,
                                        control.bgev = control.bgev),
                  control.predictor = list(compute = TRUE),
                  control.fixed = list(prec=1000),
                  control.compute = list(cpo = TRUE),
                  control.inla = list(int.strategy = "eb"),
                  verbose=FALSE, safe=TRUE)
  
  inlares2 = inla(formula2,
                  family = "bgev",
                  data = amax2,
                  control.family = list(hyper = hyper.bgev,
                                        control.bgev = control.bgev),
                  control.predictor = list(compute = TRUE),
                  control.fixed = list(prec=1000),
                  control.compute = list(cpo = TRUE),
                  control.inla = list(int.strategy = "eb"),
                  verbose=FALSE, safe=TRUE)
  
  
  eta1=rowSums(matrix(rep(t(inlares1$summary.fixed$mean),dim(amax2)[1]),nrow=dim(amax2)[1],byrow=TRUE)*as.matrix(amax2[,.(intercept,wetterdays_JJA,precip_JJA,temp_JJA,masl,lon,lat)]))
  sbeta1=inlares1$summary.hyperpar$mean[1]+inlares1$summary.hyperpar$mean[3]*amax2$lon+inlares1$summary.hyperpar$mean[4]*amax2$lat+
    inlares1$summary.hyperpar$mean[5]*amax2$temp_JJA
  xi1=inlares1$summary.hyperpar$mean[2]+inlares1$summary.hyperpar$mean[6]*amax2$lon+inlares1$summary.hyperpar$mean[7]*amax2$lat+
    inlares1$summary.hyperpar$mean[8]*amax2$temp_JJA
  
  
  eta2=rowSums(matrix(rep(t(inlares2$summary.fixed$mean),dim(amax2)[1]),nrow=dim(amax2)[1],byrow=TRUE)*as.matrix(amax2[,.(intercept,wetterdays_JJA_static_clim,precip_JJA_static_clim,temp_JJA_static_clim,masl,lon,lat)]))
  sbeta2=inlares2$summary.hyperpar$mean[1]+inlares2$summary.hyperpar$mean[3]*amax2$lon+inlares2$summary.hyperpar$mean[4]*amax2$lat+
    inlares2$summary.hyperpar$mean[5]*amax2$temp_JJA_static_clim
  xi2=inlares2$summary.hyperpar$mean[2]+inlares2$summary.hyperpar$mean[6]*amax2$lon+inlares2$summary.hyperpar$mean[7]*amax2$lat+
    inlares2$summary.hyperpar$mean[8]*amax2$temp_JJA_static_clim
  
  eta1=tail(eta1,dim(cvdata)[1])
  sbeta1=tail(sbeta1,dim(cvdata)[1])
  xi1=tail(xi1,dim(cvdata)[1])
  
  
  gevpar_estim1=giveme.gev.par(q= eta1,
                               sbeta=sbeta1,
                               alpha=p.alpha,beta=p.beta,
                               xi=xi1)
  
  eta2=tail(eta2,dim(cvdata)[1])
  sbeta2=tail(sbeta2,dim(cvdata)[1])
  xi2=tail(xi2,dim(cvdata)[1])
  
  gevpar_estim2=giveme.gev.par(q= eta2,
                               sbeta=sbeta2,
                               alpha=p.alpha,beta=p.beta,
                               xi=xi2)
  
  
  
  #-------------------------------------- Spatial model---------------------------------------------#
  formula3 = inla.mdata(y,amax2[,.(lon,lat,temp_JJA)] ,amax2[,.(lon,lat,temp_JJA)]) ~ -1 + beta0_intercept + wetterdays_JJA+precip_JJA+temp_JJA+masl+lon+lat+f(field,model=spde)
  formula4 = inla.mdata(y,amax2[,.(lon,lat,temp_JJA_static_clim)] ,amax2[,.(lon,lat,temp_JJA_static_clim)]) ~ -1 + beta0_intercept + wetterdays_JJA_static_clim+precip_JJA_static_clim+temp_JJA_static_clim + masl+lon+lat+f(field,model=spde)
  
  coords=amax2[,.(lon,lat)]
  
  max.edge=c(1,3);
  cutoff=0.25;
  offset=c(2,2)
  
  
  if(p==1){
    mesh=inla.mesh.2d(loc=coords,max.edge=max.edge,cutoff=cutoff, offset = offset)
  }
  
  Amat <- inla.spde.make.A(mesh=mesh, 
                           loc=as.matrix(coords),index=1:dim(coords)[1])
  
  spde= inla.spde2.pcmatern(mesh=mesh, alpha=2,prior.range=c(20,0.1),prior.sigma=c(2,0.1))
  s.index <- inla.spde.make.index(name="field", n.spde=spde$n.spde)
  
  stackinfo3<- inla.stack(data=list(y=amax2$y),
                          A=list(Amat,1),
                          effects=list(c(s.index,beta0_intercept=1),
                                       c(list(wetterdays_JJA=amax2$wetterdays_JJA),
                                         list(precip_JJA=amax2$precip_JJA),list(temp_JJA=amax2$temp_JJA),
                                         list(masl=amax2$masl),
                                         list(lon=amax2$lon),list(lat=amax2$lat))), tag="field")
  
  
  stackinfo4<- inla.stack(data=list(y=amax2$y),
                          A=list(Amat,1),
                          effects=list(c(s.index,beta0_intercept=1),
                                       c(list(wetterdays_JJA_static_clim=amax2$wetterdays_JJA_static_clim),
                                         list(precip_JJA_static_clim=amax2$precip_JJA_static_clim),
                                         list(temp_JJA_static_clim=amax2$temp_JJA_static_clim),
                                         list(masl=amax2$masl),
                                         list(lon=amax2$lon),list(lat=amax2$lat))), tag="field")
  
  inlares3 = inla(formula3,
                  family = "bgev",
                  data = inla.stack.data(stackinfo3),
                  control.family = list(hyper = hyper.bgev,
                                        control.bgev = control.bgev),
                  control.predictor = list(A=inla.stack.A(stackinfo3),compute=TRUE),
                  control.fixed = list(prec=1000),
                  control.compute = list(cpo = TRUE),
                  control.inla = list(int.strategy = "eb"),
                  verbose=FALSE, safe=TRUE)
  
  inlares4 = inla(formula4,
                  family = "bgev",
                  data = inla.stack.data(stackinfo4),
                  control.family = list(hyper = hyper.bgev,
                                        control.bgev = control.bgev),
                  control.predictor = list(A=inla.stack.A(stackinfo4),compute=TRUE),
                  control.fixed = list(prec=1000),
                  control.compute = list(cpo = TRUE),
                  control.inla = list(int.strategy = "eb"),
                  verbose=FALSE, safe=TRUE)
  
  #,precip_JJA,temp_JJA
  eta3=as.vector(Amat%*%inlares3$summary.random$field$mean)+inlares3$summary.fixed$mean[1]+inlares3$summary.fixed$mean[2]*amax2$wetterdays_JJA+
    inlares3$summary.fixed$mean[3]*amax2$precip_JJA+inlares3$summary.fixed$mean[4]*amax2$temp_JJA+inlares3$summary.fixed$mean[5]*amax2$masl+
    inlares3$summary.fixed$mean[6]*amax2$lon+inlares3$summary.fixed$mean[7]*amax2$lat
  sbeta3=inlares3$summary.hyperpar$mean[1]+inlares3$summary.hyperpar$mean[3]*amax2$lon+inlares3$summary.hyperpar$mean[4]*amax2$lat+
    inlares3$summary.hyperpar$mean[5]*amax2$temp_JJA
  
  xi3=inlares3$summary.hyperpar$mean[2]+inlares3$summary.hyperpar$mean[6]*amax2$lon+inlares3$summary.hyperpar$mean[7]*amax2$lat+
    inlares3$summary.hyperpar$mean[8]*amax2$temp_JJA

  
  eta4=as.vector(Amat%*%inlares4$summary.random$field$mean)+inlares4$summary.fixed$mean[1]+inlares4$summary.fixed$mean[2]*amax2$wetterdays_JJA_static_clim+
    inlares4$summary.fixed$mean[3]*amax2$precip_JJA_static_clim+inlares4$summary.fixed$mean[4]*amax2$temp_JJA_static_clim+inlares4$summary.fixed$mean[5]*amax2$masl+
    inlares4$summary.fixed$mean[6]*amax2$lon+inlares4$summary.fixed$mean[7]*amax2$lat
  
  sbeta4=inlares4$summary.hyperpar$mean[1]+inlares4$summary.hyperpar$mean[3]*amax2$lon+inlares4$summary.hyperpar$mean[4]*amax2$lat+
    inlares4$summary.hyperpar$mean[5]*amax2$temp_JJA_static_clim
  
  xi4=inlares4$summary.hyperpar$mean[2]+inlares4$summary.hyperpar$mean[6]*amax2$lon+inlares4$summary.hyperpar$mean[7]*amax2$lat+
    inlares4$summary.hyperpar$mean[8]*amax2$temp_JJA_static_clim
  
  #---------------------------------------------------------------#
  eta3=tail(eta3,dim(cvdata)[1])
  sbeta3=tail(sbeta3,dim(cvdata)[1])
  xi3=tail(xi3,dim(cvdata)[1])
  
  gevpar_estim3=giveme.gev.par(q= eta3,
                               sbeta=sbeta3,
                               alpha=p.alpha,beta=p.beta,
                               xi=xi3)
  
  eta4=tail(eta4,dim(cvdata)[1])
  sbeta4=tail(sbeta4,dim(cvdata)[1])
  xi4=tail(xi4,dim(cvdata)[1])
  
  gevpar_estim4=giveme.gev.par(q= eta4,
                               sbeta=sbeta4,
                               alpha=p.alpha,beta=p.beta,
                               xi=xi4)
  
  #---------------------------------------------------------------#
  res1=cbind(cvdata0,mu=gevpar_estim1$mu,sigma=gevpar_estim1$sigma,xi=gevpar_estim1$xi,clim_mean=y_meansd$mean_y,clim_sd=y_meansd$sd_y)
  res2=cbind(cvdata0,mu=gevpar_estim2$mu,sigma=gevpar_estim2$sigma,xi=gevpar_estim2$xi,clim_mean=y_meansd$mean_y,clim_sd=y_meansd$sd_y)
  res3=cbind(cvdata0,mu=gevpar_estim3$mu,sigma=gevpar_estim3$sigma,xi=gevpar_estim3$xi,clim_mean=y_meansd$mean_y,clim_sd=y_meansd$sd_y)
  res4=cbind(cvdata0,mu=gevpar_estim4$mu,sigma=gevpar_estim4$sigma,xi=gevpar_estim4$xi,clim_mean=y_meansd$mean_y,clim_sd=y_meansd$sd_y)
  
  
  #Frequentist model:
  par_normgev=c()
  uniquelocs=unique(amax_data[,.(lon,lat)])
  n_loc=dim(cvdata0)[1]
  for(j in 1:n_loc){
    par_normgev=rbind(par_normgev,optim(par=c(0.2,0.2,0.2),optim_gev,
                                        y=amax_data[year!=thisyear & lon==uniquelocs$lon[j] & lat==uniquelocs$lat[j],y],
                                        lower=c(-10,0.0001,0.001),
                                        upper=c(100,100,3),method="L-BFGS-B")$par)
  }
  
  colnames(par_normgev)=c("loc_gev","scale_gev","shape_gev")
  
  res1=cbind(res1,par_normgev)
  res1[,year_specific:=TRUE]
  res1[,spatial:=FALSE]
  
  inlaparam=cbind(matrix(rep(inlares1$summary.fixed$`0.5quant`,each=dim(res1)[1]),byrow=FALSE,nrow=dim(res1)[1]),
                  matrix(rep(c(inlares1$summary.hyperpar$`0.5quant`,NA,NA),each=dim(res1)[1]),byrow=FALSE,nrow=dim(res1)[1]))
  colnames(inlaparam)=c(paste0("beta0.5_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                        paste0("quant0.5_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_u=cbind(matrix(rep(inlares1$summary.fixed$`0.975quant`,each=dim(res1)[1]),byrow=FALSE,nrow=dim(res1)[1]),
                    matrix(rep(c(inlares1$summary.hyperpar$`0.975quant`,NA,NA),each=dim(res1)[1]),byrow=FALSE,nrow=dim(res1)[1]))
  colnames(inlaparam_u)=c(paste0("beta0.975_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.975_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_l=cbind(matrix(rep(inlares1$summary.fixed$`0.025quant`,each=dim(res1)[1]),byrow=FALSE,nrow=dim(res1)[1]),matrix(rep(c(inlares1$summary.hyperpar$`0.025quant`,NA,NA),each=dim(res1)[1]),byrow=FALSE,nrow=dim(res1)[1]))
  colnames(inlaparam_l)=c(paste0("beta0.025_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.025_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  res1=cbind(res1,inlaparam,inlaparam_l,inlaparam_u)
  
  #---------------------------------------------------------------#
  
  res2=cbind(res2,par_normgev)
  res2[,year_specific:=FALSE]
  res2[,spatial:=FALSE]
  
  inlaparam=cbind(matrix(rep(inlares2$summary.fixed$`0.5quant`,each=dim(res2)[1]),byrow=FALSE,nrow=dim(res1)[1]),matrix(rep(c(inlares2$summary.hyperpar$`0.5quant`,NA,NA),each=dim(res2)[1]),byrow=FALSE,nrow=dim(res2)[1]))
  colnames(inlaparam)=c(paste0("beta0.5_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                        paste0("quant0.5_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_u=cbind(matrix(rep(inlares2$summary.fixed$`0.975quant`,each=dim(res2)[1]),byrow=FALSE,nrow=dim(res2)[1]),matrix(rep(c(inlares2$summary.hyperpar$`0.975quant`,NA,NA),each=dim(res2)[1]),byrow=FALSE,nrow=dim(res2)[1]))
  colnames(inlaparam_u)=c(paste0("beta0.975_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.975_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_l=cbind(matrix(rep(inlares2$summary.fixed$`0.025quant`,each=dim(res2)[1]),byrow=FALSE,nrow=dim(res2)[1]),matrix(rep(c(inlares2$summary.hyperpar$`0.025quant`,NA,NA),each=dim(res2)[1]),byrow=FALSE,nrow=dim(res2)[1]))
  colnames(inlaparam_l)=c(paste0("beta0.025_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.025_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  
  res2=cbind(res2,inlaparam,inlaparam_l,inlaparam_u)
  
  #---------------------------------------------------------------#
  
  
  res3=cbind(res3,par_normgev)
  res3[,year_specific:=TRUE]
  res3[,spatial:=TRUE]
  
  inlaparam=cbind(matrix(rep(inlares3$summary.fixed$`0.5quant`,each=dim(res3)[1]),byrow=FALSE,nrow=dim(res3)[1]),matrix(rep(c(inlares3$summary.hyperpar$`0.5quant`),each=dim(res3)[1]),byrow=FALSE,nrow=dim(res3)[1]))
  colnames(inlaparam)=c(paste0("beta0.5_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                        paste0("quant0.5_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_u=cbind(matrix(rep(inlares3$summary.fixed$`0.975quant`,each=dim(res3)[1]),byrow=FALSE,nrow=dim(res3)[1]),matrix(rep(c(inlares3$summary.hyperpar$`0.975quant`),each=dim(res3)[1]),byrow=FALSE,nrow=dim(res3)[1]))
  colnames(inlaparam_u)=c(paste0("beta0.975_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.975_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_l=cbind(matrix(rep(inlares3$summary.fixed$`0.025quant`,each=dim(res3)[1]),byrow=FALSE,nrow=dim(res3)[1]),matrix(rep(c(inlares3$summary.hyperpar$`0.025quant`),each=dim(res3)[1]),byrow=FALSE,nrow=dim(res3)[1]))
  colnames(inlaparam_l)=c(paste0("beta0.025_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.025_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  
  res3=cbind(res3,inlaparam,inlaparam_l,inlaparam_u)
  
  
  
  res4=cbind(res4,par_normgev)
  res4[,year_specific:=FALSE]
  res4[,spatial:=TRUE]
  inlaparam=cbind(matrix(rep(inlares4$summary.fixed$`0.5quant`,each=dim(res4)[1]),byrow=FALSE,nrow=dim(res3)[1]),matrix(rep(c(inlares4$summary.hyperpar$`0.5quant`),each=dim(res4)[1]),byrow=FALSE,nrow=dim(res4)[1]))
  colnames(inlaparam)=c(paste0("beta0.5_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                        paste0("quant0.5_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_u=cbind(matrix(rep(inlares4$summary.fixed$`0.975quant`,each=dim(res4)[1]),byrow=FALSE,nrow=dim(res4)[1]),matrix(rep(c(inlares4$summary.hyperpar$`0.975quant`),each=dim(res4)[1]),byrow=FALSE,nrow=dim(res4)[1]))
  colnames(inlaparam_u)=c(paste0("beta0.975_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.975_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  inlaparam_l=cbind(matrix(rep(inlares4$summary.fixed$`0.025quant`,each=dim(res4)[1]),byrow=FALSE,nrow=dim(res4)[1]),matrix(rep(c(inlares4$summary.hyperpar$`0.025quant`),each=dim(res4)[1]),byrow=FALSE,nrow=dim(res4)[1]))
  colnames(inlaparam_l)=c(paste0("beta0.025_",c("intercept","wetterdays_JJA","precip_JJA","temp_JJA","masl","lon","lat")),
                          paste0("quant0.025_",c("spread","tail","beta1_spread","beta2_spread","beta3_spread","beta4_tail","beta5_tail","beta6_tail","range","sig")))
  
  
  res4=cbind(res4,inlaparam,inlaparam_l,inlaparam_u)
  
  reslist[[p]]=rbind(res1,res2,res3,res4)
  p=p+1
}

reslist=rbindlist(reslist)

fwrite(reslist,"/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/bGev_tutorial/inlares_advanced.csv")
