library(INLA)
library(evgam)
library(evd)
library(data.table)
loc_thea = "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/"
setwd(loc_thea)
#INLA:::inla.binary.install()
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

#try with precip data:
loc_thea = "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/"
setwd(loc_thea)

#Upload hourly data (change 60 to 1440 to get daily data):
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")[,.(year,masl,stid,lon,lat,
                                                               wetterdays_JJA,wetterdays_annual,
                                                               precip_JJA,precip_annual,
                                                               temp_JJA,temp_annual,y)]


amax2= data.table(scale(amax_data[,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl)]))
amax2[,intercept:=1]

std_data=amax_data[,.(mean_y=mean(y),sd_y=sd(y)),]

#---------------- With and without year-specific covariates -----------------------------#
lon=amax2$lon;lat=amax2$lat
formula1 = inla.mdata(y,cbind(lon,lat) ,cbind(lon,lat)) ~ -1 + intercept + wetterdays_JJA+precip_JJA+temp_JJA+masl+lon+lat
formula2 = inla.mdata(y,cbind(lon,lat),cbind(lon,lat)) ~ -1 + intercept +masl+lon+lat

inlares1 = inla(formula1,
                family = "bgev",
                data = amax2,
                control.family = list(hyper = hyper.bgev,
                                      control.bgev = control.bgev),
                control.predictor = list(compute = TRUE),
                control.fixed = list(prec=1000),
                control.compute = list(cpo = TRUE),
                #control.inla = list(int.strategy = "eb"),
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
sbeta1=inlares1$summary.hyperpar$mean[1]+inlares1$summary.hyperpar$mean[3]*amax2$lon+inlares1$summary.hyperpar$mean[4]*amax2$lat
xi1=inlares1$summary.hyperpar$mean[2]+inlares1$summary.hyperpar$mean[5]*amax2$lon+inlares1$summary.hyperpar$mean[6]*amax2$lat

eta2=rowSums(matrix(rep(t(inlares2$summary.fixed$mean),dim(amax2)[1]),nrow=dim(amax2)[1],byrow=TRUE)*as.matrix(amax2[,.(intercept,masl,lon,lat)]))
sbeta2=inlares2$summary.hyperpar$mean[1]+inlares2$summary.hyperpar$mean[3]*amax2$lon+inlares2$summary.hyperpar$mean[4]*amax2$lat
xi2=inlares2$summary.hyperpar$mean[2]+inlares2$summary.hyperpar$mean[5]*amax2$lon+inlares2$summary.hyperpar$mean[6]*amax2$lat


gevpar_estim1=giveme.gev.par(q= eta1,
                             sbeta=sbeta1,
                             alpha=p.alpha,beta=p.beta,
                             xi=xi1)


gevpar_estim2=giveme.gev.par(q= eta2,
                             sbeta=sbeta2,
                             alpha=p.alpha,beta=p.beta,
                             xi=xi2)



#-------------------------------------- Spatial model---------------------------------------------#
formula3 = inla.mdata(y,amax2[,.(lon,lat)] ,amax2[,.(lon,lat)]) ~ -1 + beta0_intercept + wetterdays_JJA+precip_JJA+temp_JJA+masl+lon+lat+f(field,model=spde)
formula4 = inla.mdata(y,amax2[,.(lon,lat)] ,amax2[,.(lon,lat)]) ~ -1 + beta0_intercept + masl+lon+lat+f(field,model=spde)


coords=amax2[,.(lon,lat)]
max.edge=c(1,3);
cutoff=0.05;
offset=c(2,2)
mesh=inla.mesh.2d(loc=coords,max.edge=max.edge,cutoff=cutoff, offset = offset)
plot(mesh)
points(coords)


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
                                     c(list(masl=amax2$masl),
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


eta4=as.vector(Amat%*%inlares4$summary.random$field$mean)+inlares4$summary.fixed$mean[1]+inlares4$summary.fixed$mean[2]*amax2$masl+
  inlares4$summary.fixed$mean[3]*amax2$lon+inlares4$summary.fixed$mean[4]*amax2$lat
sbeta4=inlares4$summary.hyperpar$mean[1]+inlares4$summary.hyperpar$mean[3]*amax2$lon+inlares4$summary.hyperpar$mean[4]*amax2$lat
xi4=inlares4$summary.hyperpar$mean[2]+inlares4$summary.hyperpar$mean[5]*amax2$lon+inlares4$summary.hyperpar$mean[6]*amax2$lat

gevpar_estim4=giveme.gev.par(q= eta4,
                             sbeta=sbeta4,
                             alpha=p.alpha,beta=p.beta,
                             xi=xi4)

#----------------------------------------------------------------------------------------------------#
#Frequentist model:
par_normgev=c()
uniquelocs=unique(amax2[,.(lon,lat)])
n_loc=dim(uniquelocs)[1]
for(j in 1:n_loc){
  par_normgev=rbind(par_normgev,optim(par=c(0.2,0.2,0.2),optim_gev,
                                      y=amax2[lon==uniquelocs$lon[j] & lat==uniquelocs$lat[j],y],
                                      lower=c(-10,0.0001,0.0001),
                                      upper=c(10,100,1),method="L-BFGS-B")$par)
}
#----------------------------------------------------------------------------------------------------#



#Location specifics
par(mfrow=c(2,2))
uniquelocs=unique(amax_data[,.(lon,lat)])
resmat_gev=c()
resmat_bgev=c()
resmat_spatgev=c()
for(locnum in 1:n_loc){
  print(locnum)
  # Plot results
  return_periods=c(2,5,10,25,50)
  quantiles=1-1/return_periods
  thisdata=amax_data[lon==uniquelocs$lon[locnum] & lat==uniquelocs$lat[locnum],]
  thisdata[,gposition:=lapply(.SD,function(x) 1/(1-gringorten(x))), .SDcols = "y"]
  
  bgevres=c()
  gevres=c()
  meaninfo=std_data
  spatgevres=c()
  
  for(j in 1:length(quantiles)){
    bgevres[j]=qbgev(quantiles[j],mu=unique(gevpar_estim2$mu)[locnum],
                     sigma=unique(gevpar_estim2$sigma)[locnum],
                     xi=unique(gevpar_estim2$xi)[locnum],p_b=0.2,s=5)
    
    spatgevres[j]=qbgev(quantiles[j],mu=unique(gevpar_estim4$mu)[locnum],
                        sigma=unique(gevpar_estim4$sigma)[locnum],
                        xi=unique(gevpar_estim4$xi)[locnum],p_b=0.2,s=5)
    
    gevres[j]=evd::qgev(quantiles[j],par_normgev[locnum,1],par_normgev[locnum,2],par_normgev[locnum,3])
  }
  
  plot(return_periods,bgevres*meaninfo$sd_y+meaninfo$mean_y,type="o",ylim=c(0,40),xlab="Return period",ylab="mm");grid()
  lines(return_periods,return_level_gev(return_periods,mu=unique(gevpar_estim2$mu)[locnum],
                                        sigma=unique(gevpar_estim2$sigma)[locnum],
                                        xi=unique(gevpar_estim2$xi)[locnum])*meaninfo$sd_y+meaninfo$mean_y,col="red")
  lines(return_periods,return_level_bgev(return_periods,mu=unique(gevpar_estim2$mu)[locnum],
                                         sigma=unique(gevpar_estim2$sigma)[locnum],
                                         xi=unique(gevpar_estim2$xi)[locnum],
                                         p_a=p.alpha,p_b=p.beta,s=5)*meaninfo$sd_y+meaninfo$mean_y,col="yellow",lty=2)
  lines(return_periods,return_level_bgev(return_periods,mu=unique(gevpar_estim4$mu)[locnum],
                                         sigma=unique(gevpar_estim4$sigma)[locnum],
                                         xi=unique(gevpar_estim4$xi)[locnum],
                                         p_a=p.alpha,p_b=p.beta,s=5)*meaninfo$sd_y+meaninfo$mean_y,col="orange",type="o",lty=2)
  
  lines(return_periods,gevres*meaninfo$sd_y+meaninfo$mean_y,col="blue",lty=2,type="o")
  points(thisdata$gposition,thisdata$y,pch=3)
  title(paste0("Lon ",round(uniquelocs$lon[locnum],1),", lat ", round(uniquelocs$lat[locnum],1)," and nyears ",dim(thisdata)[1]))
  
  
  resmat_bgev=rbind(resmat_bgev,return_level_bgev(return_periods,mu=unique(gevpar_estim2$mu)[locnum],
                                                  sigma=unique(gevpar_estim2$sigma)[locnum],
                                                  xi=unique(gevpar_estim2$xi)[locnum],
                                                  p_a=p.alpha,p_b=p.beta,s=5)*meaninfo$sd_y+meaninfo$mean_y)
  
  resmat_spatgev=rbind(resmat_spatgev,return_level_bgev(return_periods,mu=unique(gevpar_estim4$mu)[locnum],
                                                        sigma=unique(gevpar_estim4$sigma)[locnum],
                                                        xi=unique(gevpar_estim4$xi)[locnum],
                                                        p_a=p.alpha,p_b=p.beta,s=5)*meaninfo$sd_y+meaninfo$mean_y)
  
  
  resmat_gev=rbind(resmat_gev,gevres*meaninfo$sd_y+meaninfo$mean_y)
}



par(mfrow=c(1,2))
plot(resmat_gev,resmat_bgev); lines(c(-100,100),c(-100,100),col="red",xlim=c(0,140),ylim=c(0,140));grid()
plot(resmat_spatgev,resmat_bgev); lines(c(-100,100),c(-100,100),col="red");grid()

field=inlares3$summary.random$field$mean
proj = inla.mesh.projector(mesh,dims=c(300, 300))
field.proj = inla.mesh.project(proj,field)
image.plot(list(x=proj$x,y=proj$y,z=field.proj))
points(coords,col="gray",cex=0.2)


field=inlares4$summary.random$field$mean
proj = inla.mesh.projector(mesh,dims=c(300, 300))
field.proj = inla.mesh.project(proj,field)
image.plot(list(x=proj$x,y=proj$y,z=field.proj))
points(coords,col="gray",cex=0.2)

