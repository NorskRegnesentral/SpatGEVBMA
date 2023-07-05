library(INLA)
library(evgam)
library(evd)
library(data.table)
loc_thea = "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/"
setwd(loc_thea)
#INLA:::inla.binary.install()

giveme.gev.par = function(q, sbeta, alpha, beta, xi)
{
  .mu = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    c = -log(alpha)
    if (all(xi > 0.0)) {
      tmp0 = (c^(-xi) - 1)/xi
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(q - (sbeta/dbeta) * tmp0)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      tmp0 = log(c)
      return(q + (sbeta/dbeta) * tmp0)
    } else {
      stop("mixed case not implemented")
    }
  }
  .sigma = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    if (all(xi > 0.0)) {
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(sbeta/dbeta)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      return(sbeta/dbeta)
    } else {
      stop("mixed case not implemented")
    }
  }
  return(list(mu = .mu(q, sbeta, alpha, beta, xi),
              sigma = .sigma(q, sbeta, alpha, beta, xi),
              xi = xi))
}


map.tail = function(x, interval, inverse = FALSE) {
  if (!inverse) {
    return (interval[1] + (interval[2] - interval[1]) * exp(x)/(1.0 + exp(x)))
  } else {
    return (log((x-interval[1])/(interval[2]-x)))
  }
}

#---------------------------------------------------------------------------------------------------------------------------#

# Example 1: Spread and tail not dependent on covariates.
set.seed(2)
n_loc=4
n_each=1000
n_tot=n_loc*n_each
x = c(rep(rnorm(n_loc,sd=0.5),each=n_each))
eta.x = 1 + 0.4*x

spread = 0.3
tail = 0.1

p.alpha = 0.5
p.beta = 0.25

par = giveme.gev.par(q=eta.x,sbeta=spread,alpha=p.alpha,beta=p.beta,xi=tail)

y=numeric(n_tot)
for(i in 1:n_tot){
  y[i]=evd::rgev(1,loc=par$mu[i],scale=par$sigma,shape=par$xi)
}

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
                    q.mix= c(0.05, 0.20),
                    # the Beta(s1, s2) mixing distribution parameters.
                    # Hard-coded for the moment: s1=s2=5
                    beta.ab = 5)


# INLA fit:
null.matrix = matrix(nrow = n_tot, ncol= 0) #empty
spread.x = null.matrix #indicate that spread and tail are covariate free.
tail.x = null.matrix

data.bgev = data.frame(y = y, intercept = 1, x = x, spread.x = spread.x, tail.x = tail.x)
formula = inla.mdata(y, spread.x, tail.x) ~ -1 + intercept + x

r1 = inla(formula,
          family = "bgev",
          data = data.bgev,
          control.family = list(hyper = hyper.bgev,
                                control.bgev = control.bgev),
          control.predictor = list(compute = TRUE),
          control.fixed = list(prec=1000),
          control.compute = list(cpo = TRUE),
          control.inla = list(int.strategy = "eb"),
          verbose=FALSE, safe=TRUE)

round(r1$summary.fixed,4)
round(r1$summary.hyperpar,4)

gevpar_estim=giveme.gev.par(q=r1$summary.fixed$mean[1]+r1$summary.fixed$mean[2]*x,
               sbeta=r1$summary.hyperpar$mean[1],
               alpha=p.alpha,beta=p.beta,xi=r1$summary.hyperpar$mean[2])

par_normgev=c()
for(j in 1:n_loc){
  par_normgev=rbind(par_normgev,optim(par=c(1,1,1),optim_gev,
                  y=y[1:n_each+n_each*(j-1)],lower=c(0.0001,0.0001,-5),
                  upper=c(100,100,5),method="L-BFGS-B")$par)
}

#-------------------------------------------------------------------------------------------#
par(mfrow=c(2,2))
for(locnum in 1:n_loc){
  # Plot results

  return_periods=c(2,5,10,25,50,100,200)
  quantiles=1-1/return_periods
  
  bgevres=c()
  gevres=c()
  for(j in 1:length(quantiles)){
    bgevres[j]=qbgev(quantiles[j],mu=gevpar_estim$mu[(locnum-1)*n_each+1],
                    sigma=gevpar_estim$sigma,xi=gevpar_estim$xi,p_b=0.2,s=5)
    gevres[j]=evd::qgev(quantiles[j],par_normgev[locnum,1],par_normgev[locnum,2],par_normgev[locnum,3])
  }
  
  plot(return_periods,bgevres,type="o",ylim=c(1,2));grid()
  lines(return_periods,return_level_gev(return_periods,mu=gevpar_estim$mu[(locnum-1)*n_each+1],
                                        sigma=gevpar_estim$sigma,xi=gevpar_estim$xi),col="red")
  lines(return_periods,gevres,col="blue",lty=2,type="o")
}



#-------------------------------------------------------------------------------------------#
# Example 2:
p.alpha = 0.5
p.beta = 0.25

n_loc=4
n_each=1000
n_tot=n_loc*n_each

x1 = c(rep(rnorm(n_loc,mean=10,sd=2)/10,each=n_each))

eta.x = 1 + 0.4*x1

x2 = c(rep(rnorm(n_loc, sd= 1,mean=3),each=n_each))
s.x = exp(0.1 + 0.3*x2)

x3 = rep(runif(n_loc,-0.25,2),each=n_each)

t.x = 0.1 + 0.2*x3
tail.intern = map.tail(t.x, tail.interval, inverse=TRUE)

par = giveme.gev.par(q = eta.x, sbeta = s.x, alpha = p.alpha, beta = p.beta,
                     xi = t.x)

y = numeric(n_tot)
for(i in 1:n_tot)
  y[i] = evd::rgev(1, loc = par$mu[i], scale = par$sigma[i], shape = par$xi[i])


hyper.beta1 = hyper.beta2 = list(prior = "normal",
                                 param = c(0, 300),
                                 initial = 0)


hyper.bgev = list(spread = hyper.spread,
                  tail = hyper.tail,
                  beta1 = hyper.beta1,
                  beta2 = hyper.beta2)

spread.x = x2
tail.x = x3
formula = inla.mdata(y, spread.x, tail.x) ~ -1 + intercept + x1
data.bgev = data.frame(y = y, intercept = 1, x1 = x1, spread.x = spread.x, tail.x = tail.x)

r2 = inla(formula,
          family = "bgev",
          data = data.bgev,
          control.family = list(hyper = hyper.bgev,
                                control.bgev = control.bgev),
          control.predictor = list(compute = TRUE),
          control.fixed = list(prec=1000),
          control.compute = list(cpo = TRUE),
          control.inla = list(int.strategy = "eb"),
          verbose=FALSE, safe=TRUE)

round(r2$summary.fixed,4)
round(r2$summary.hyperpar,4)


gevpar_estim=giveme.gev.par(q=r2$summary.fixed$mean[1]+r2$summary.fixed$mean[2]*x1,
                            sbeta=r2$summary.hyperpar$mean[1]+r2$summary.hyperpar$mean[3]*x2,
                            alpha=p.alpha,beta=p.beta,xi=r2$summary.hyperpar$mean[2]+
                          r2$summary.hyperpar$mean[4]*x3)

par_normgev=c()
for(j in 1:n_loc){
  par_normgev=rbind(par_normgev,optim(par=c(1,1,1),optim_gev,
                                      y=y[1:n_each+n_each*(j-1)],
                                      lower=c(0.0001,0.0001,-5),
                                      upper=c(100,100,5),method="L-BFGS-B")$par)
}


par(mfrow=c(2,2))
for(locnum in 1:n_loc){
  # Plot results
  return_periods=c(2,5,10,25,50,100)
  quantiles=1-1/return_periods
  
  bgevres=c()
  gevres=c()
  for(j in 1:length(quantiles)){
    bgevres[j]=qbgev(quantiles[j],mu=gevpar_estim$mu[(locnum-1)*n_each+1],
                     sigma=gevpar_estim$sigma[(locnum-1)*n_each+1],
                     xi=gevpar_estim$xi[(locnum-1)*n_each+1],p_b=0.2,s=5)
    gevres[j]=evd::qgev(quantiles[j],par_normgev[locnum,1],par_normgev[locnum,2],par_normgev[locnum,3])
  }
  
  plot(return_periods,bgevres,type="o",ylim=c(0,10));grid()
  lines(return_periods,return_level_gev(return_periods,mu=gevpar_estim$mu[(locnum-1)*n_each+1],
                                        sigma=gevpar_estim$sigma[(locnum-1)*n_each+1],xi=gevpar_estim$xi[(locnum-1)*n_each+1]),col="red")
  lines(return_periods,return_level_bgev(return_periods,mu=gevpar_estim$mu[(locnum-1)*n_each+1],
                                         sigma=gevpar_estim$sigma[(locnum-1)*n_each+1],xi=gevpar_estim$xi[(locnum-1)*n_each+1],
                                         p_a=p.alpha,p_b=p.beta,s=5),col="yellow",lty=2)
  lines(return_periods,gevres,col="blue",lty=2,type="o")
}

#------------------------------------------------------------------------------------------------#
#try with precip data:
loc_thea = "/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/"
setwd(loc_thea)

#Upload hourly data (change 60 to 1440 to get daily data):
amax_data=fread(file="scripts/dev_I4C/Data/AM1440_cov.csv")[,.(year,masl,stid,lon,lat,
                                                             wetterdays_JJA,wetterdays_annual,
                                                             precip_JJA,precip_annual,
                                                             temp_JJA,temp_annual,y)]

amax2= data.table(scale(amax_data[,.(wetterdays_JJA,precip_JJA,y,lon,lat,temp_JJA,masl)]))
amax2[,intercept:=1]

std_data=amax_data[,.(mean_y=mean(y),sd_y=sd(y)),.(lon,lat)]


#---------------- With and without year-specific covariates -----------------------------#
formula1 = inla.mdata(y,amax2[,.(lon,lat)] ,amax2[,.(lon,lat)]) ~ -1 + intercept + wetterdays_JJA+precip_JJA+temp_JJA+masl+lon+lat
formula2 = inla.mdata(y,amax2[,.(lon,lat)] ,amax2[,.(lon,lat)]) ~ -1 + intercept +masl+lon+lat

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

#--------------------------------------------------------------------#
eta1=as.vector(matrix(rep(t(inlares1$summary.fixed$mean),dim(amax2)[1]),nrow=dim(amax2)[1],byrow=TRUE)*as.matrix(amax2[,.(intercept,wetterdays_JJA,precip_JJA,temp_JJA,masl,lon,lat)]))
sbeta1=inlares1$summary.hyperpar$mean[1]+inlares1$summary.hyperpar$mean[3]*amax2$lon+inlares1$summary.hyperpar$mean[4]*amax2$lat
xi1=inlares1$summary.hyperpar$mean[2]+inlares1$summary.hyperpar$mean[5]*amax2$lon+inlares1$summary.hyperpar$mean[6]*amax2$lat

gevpar_estim1=giveme.gev.par(q= eta1,
                            sbeta=sbeta1,
                            alpha=p.alpha,beta=p.beta,
                            xi=xi1)

eta2=as.vector(matrix(rep(t(inlares2$summary.fixed$mean),dim(amax2)[1]),nrow=dim(amax2)[1],byrow=TRUE)*as.matrix(amax2[,.(intercept,masl,lon,lat)]))
sbeta2=inlares2$summary.hyperpar$mean[1]+inlares2$summary.hyperpar$mean[3]*amax2$lon+inlares2$summary.hyperpar$mean[4]*amax2$lat
xi2=inlares2$summary.hyperpar$mean[2]+inlares2$summary.hyperpar$mean[5]*amax2$lon+inlares2$summary.hyperpar$mean[6]*amax2$lat

gevpar_estim2=giveme.gev.par(q= eta2,
                             sbeta=sbeta2,
                             alpha=p.alpha,beta=p.beta,
                             xi=xi2)


#-------------------------------------- Spatial model---------------------------------------------#
formula1 = inla.mdata(y,amax2[,.(lon,lat)] ,amax2[,.(lon,lat)]) ~ -1 + beta0_intercept + wetterdays_JJA+precip_JJA+temp_JJA+masl+lon+lat+f(field,model=spde)
formula2 = inla.mdata(y,amax2[,.(lon,lat)] ,amax2[,.(lon,lat)]) ~ -1 + beta0_intercept + masl+lon+lat+f(field,model=spde)

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

stackinfo1<- inla.stack(data=list(y=amax2$y),
                        A=list(Amat,1),
                        effects=list(c(s.index,beta0_intercept=1),
                                     c(list(wetterdays_JJA=amax2$wetterdays_JJA),
                                       list(precip_JJA=amax2$precip_JJA),list(temp_JJA=amax2$temp_JJA),
                                       list(masl=amax2$masl),
                                       list(lon=amax2$lon),list(lat=amax2$lat))), tag="field")


stackinfo2<- inla.stack(data=list(y=amax2$y),
                        A=list(Amat,1),
                        effects=list(c(s.index,beta0_intercept=1),
                                     c(list(masl=amax2$masl),
                                       list(lon=amax2$lon),list(lat=amax2$lat))), tag="field")


inlares3 = inla(formula1,
                family = "bgev",
                data = inla.stack.data(stackinfo1),
                control.family = list(hyper = hyper.bgev,
                                      control.bgev = control.bgev),
                control.predictor = list(A=inla.stack.A(stackinfo1),compute=TRUE),
                control.fixed = list(prec=1000),
                control.compute = list(cpo = TRUE),
                control.inla = list(int.strategy = "eb"),
                verbose=FALSE, safe=TRUE)

inlares4 = inla(formula2,
                family = "bgev",
                data = inla.stack.data(stackinfo2),
                control.family = list(hyper = hyper.bgev,
                                      control.bgev = control.bgev),
                control.predictor = list(A=inla.stack.A(stackinfo2),compute=TRUE),
                control.fixed = list(prec=1000),
                control.compute = list(cpo = TRUE),
                control.inla = list(int.strategy = "eb"),
                verbose=FALSE, safe=TRUE)


#inlares3$summary.fixed
#inlares3$summary.hyperpar

inlares4$summary.fixed
inlares4$summary.hyperpar

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

#Location specifics
par(mfrow=c(2,2))
uniquelocs=unique(std_data[,.(lon,lat)])
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
  meaninfo=std_data[lon==uniquelocs$lon[locnum] & lat==uniquelocs$lat[locnum]]
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
  
  plot(return_periods,bgevres*meaninfo$sd_y+meaninfo$mean_y,type="o",ylim=c(0,35),xlab="Return period",ylab="mm");grid()
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


#----------------------------------------------------------------------------------------------------------------------------#

