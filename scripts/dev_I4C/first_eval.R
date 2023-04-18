library(SpatGEVBMA)
library(MASS)

makeLinPred=function(R,burnin=1e2,X,include_tau=FALSE){
  #kappa could be logged.
  
  n=dim(R$TAU)[1]
  
  mu=list(); kappa=list(); xi=list()
  for(stedind in 1:dim(X)[1]){
    #Param nr 1 (mu?):
    mu[[stedind]]=rowSums((R$THETA[burnin:n,,1])*as.double(matrix(rep(X[stedind,],length(burnin:n)),byrow=TRUE,ncol=length(X[stedind,]))))+
      R$TAU[burnin:n,stedind,1]*include_tau
    
    kappa[[stedind]]=rowSums((R$THETA[burnin:n,,2])*as.double(matrix(rep(X[stedind,],length(burnin:n)),byrow=TRUE,ncol=length(X[stedind,]))))+
      R$TAU[burnin:n,stedind,2]*include_tau
    
    
    xi[[stedind]]=rowSums((R$THETA[burnin:n,,3])*as.double(matrix(rep(X[stedind,],length(burnin:n)),byrow=TRUE,ncol=length(X[stedind,]))))+
      R$TAU[burnin:n,stedind,3]*include_tau
  }
  
  if(R$log.kappa==TRUE){
    kappa=lapply(kappa,exp)
  }
  
  return(list(kappa=kappa,mu=mu,xi=xi))
}

gevwrap=function(paramvec,retper=c(2,5,10,25,50,100,200)){
  prob=1-1/retper
  return(qgev(prob,loc=paramvec[1],scale=1/paramvec[2],shape=paramvec[3],lower.tail=TRUE))
}

optim_gev=function(paramvec,y){
  loc=paramvec[1]
  scale=paramvec[2]
  shape=paramvec[3]
  
  liksum=0
  for(j in 1:length(y)){
    liksum=liksum+log(dgev(y[j],loc=loc,scale=scale,shape=shape)+0.001)
  }
  
  return(-liksum)
}

#----------Frequentist simple GEV-fit---------------------#
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")[lon<7&lat<62 & lat>58 & year>=1972 & year<=2015] 
uniquelocs2=unique(amax_data[,.(lon,lat)])

freqres=c()
for(j in 1:dim(uniquelocs2)[1]){
  sub=amax_data[lon==uniquelocs2$lon[j] & lat==uniquelocs2$lat[j],.(stid,year,lon,lat,y)]
  
  freqres=rbind(freqres,optim(par=c(1,1,1),optim_gev,y=sub$y,lower=c(5,0.0001,-2),upper=c(50,50,2),method="L-BFGS-B")$par)
  
}
freqres[,2]=1/freqres[,2]

#----------------------------------------------------------------#


summarize_results=function(mcmc_res,amax_data,fit_type="annual",include_tau=1,clim_res=FALSE,burnin=10000){
  #fit_type="annual" or "annual_clim" or "full spatial"
  if(fit_type=="annual"){
    X=copy(amax_data[lon<7&lat<62 & lat>58 & year>=1972 & year<=2015,.(year,masl,stid,lon,lat,
                           wetterdays_JJA,wetterdays_annual,
                           precip_JJA,precip_annual,
                           temp_JJA,temp_annual,y)])
    
  }
  
  if(fit_type=="full_spatial"){
    X=unique(copy(amax_data[lon<7&lat<62 & lat>58 & year>=1972 & year<=2015,.(masl,lon,lat,
                                                                              wetterdays_JJA_static_clim,wetterdays_annual_static_clim,
                                                                              precip_JJA_static_clim,precip_annual_static_clim,
                                                                              temp_JJA_static_clim,temp_annual_static_clim)]))
    spatgev_data=list()
    X_unscaled=X
    X=scale(X)
    
    X=cbind(1,X)
    spatgev_data$X=X
  }
  
  if(fit_type=="annual_clim"){
    X=copy(amax_data[lon<7&lat<62 & lat>58 & year>=1972 & year<=2015,.(year,masl,stid,lon,lat,
                                                                       wetterdays_JJA_clim,wetterdays_annual_clim,
                                                                       precip_JJA_clim,precip_annual_clim,
                                                                       temp_JJA_clim,temp_annual_clim,y)])
  }
  
  if(fit_type!="full_spatial"){
    X_unscaled=copy(X)
    yinfo=X_unscaled[,.(year,stid,lon,lat,y)]
    X_unscaled[,year:=NULL];X_unscaled[,stid:=NULL];X_unscaled[,y:=NULL]
    
    spatgev_data=make_temporal_spatgev_data(X,TRUE)
    X=spatgev_data$X
  }
  
  
  if(clim_res==TRUE){
    X_means=c(colMeans(X_unscaled))
    X_sd=c(apply(X_unscaled,2,sd))
    X_clim=unique(copy(amax_data[lon<7&lat<62 & lat>58 & year>=1972 & year<=2015,.(masl,lon,lat,
                                                                                  wetterdays_JJA_static_clim,wetterdays_annual_static_clim,
                                                                                  precip_JJA_static_clim,precip_annual_static_clim,
                                                                                  temp_JJA_static_clim,temp_annual_static_clim)]))
    
    X_clim_unscaled=copy(X_clim)
    
    if(fit_type=="annual"){
      X_clim[,masl:=(masl-X_means[which(names(X_means)=="masl")])/X_sd[which(names(X_means)=="masl")]]
      X_clim[,lon:=(lon-X_means[which(names(X_means)=="lon")])/X_sd[which(names(X_means)=="lon")]]  
      X_clim[,lat:=(lat-X_means[which(names(X_means)=="lat")])/X_sd[which(names(X_means)=="lat")]]
      X_clim[,wetterdays_JJA_static_clim:=(wetterdays_JJA_static_clim-X_means[which(names(X_means)=="wetterdays_JJA")])/X_sd[which(names(X_means)=="wetterdays_JJA")]]
      X_clim[,wetterdays_annual_static_clim:=(wetterdays_annual_static_clim-X_means[which(names(X_means)=="wetterdays_annual")])/X_sd[which(names(X_means)=="wetterdays_annual")]]
      X_clim[,precip_JJA_static_clim:=(precip_JJA_static_clim-X_means[which(names(X_means)=="precip_JJA")])/X_sd[which(names(X_means)=="precip_JJA")]]
      X_clim[,precip_annual_static_clim:=(precip_annual_static_clim-X_means[which(names(X_means)=="precip_annual")])/X_sd[which(names(X_means)=="precip_annual")]]
      X_clim[,temp_JJA_static_clim:=(temp_JJA_static_clim-X_means[which(names(X_means)=="temp_JJA")])/X_sd[which(names(X_means)=="temp_JJA")]]
      X_clim[,temp_annual_static_clim:=(temp_annual_static_clim-X_means[which(names(X_means)=="temp_annual")])/X_sd[which(names(X_means)=="temp_annual")]]
    }
    
    if(fit_type=="annual_clim"){
      X_clim[,masl:=(masl-X_means[which(names(X_means)=="masl")])/X_sd[which(names(X_means)=="masl")]]
      X_clim[,lon:=(lon-X_means[which(names(X_means)=="lon")])/X_sd[which(names(X_means)=="lon")]]  
      X_clim[,lat:=(lat-X_means[which(names(X_means)=="lat")])/X_sd[which(names(X_means)=="lat")]]
      X_clim[,wetterdays_JJA_static_clim:=(wetterdays_JJA_static_clim-X_means[which(names(X_means)=="wetterdays_JJA_clim")])/X_sd[which(names(X_means)=="wetterdays_JJA_clim")]]
      X_clim[,wetterdays_annual_static_clim:=(wetterdays_annual_static_clim-X_means[which(names(X_means)=="wetterdays_annual_clim")])/X_sd[which(names(X_means)=="wetterdays_annual_clim")]]
      X_clim[,precip_JJA_static_clim:=(precip_JJA_static_clim-X_means[which(names(X_means)=="precip_JJA_clim")])/X_sd[which(names(X_means)=="precip_JJA_clim")]]
      X_clim[,precip_annual_static_clim:=(precip_annual_static_clim-X_means[which(names(X_means)=="precip_annual_clim")])/X_sd[which(names(X_means)=="precip_annual_clim")]]
      X_clim[,temp_JJA_static_clim:=(temp_JJA_static_clim-X_means[which(names(X_means)=="temp_JJA_clim")])/X_sd[which(names(X_means)=="temp_JJA_clim")]]
      X_clim[,temp_annual_static_clim:=(temp_annual_static_clim-X_means[which(names(X_means)=="temp_annual_clim")])/X_sd[which(names(X_means)=="temp_annual_clim")]]
      
    }

    X_clim=cbind(1,X_clim) 

    
    
  }else{
    X_clim=NULL
    X_clim_unscaled=NULL
    
    if(fit_type=="full_spatial"){
      
      X_clim=X
      X_clim_unscaled=X_unscaled
    }
  }
  
  if(clim_res==TRUE){
    linpred=makeLinPred(mcmc_res,burnin=burnin,X_clim,include_tau=include_tau)

  }else{
    linpred=makeLinPred(mcmc_res,burnin=burnin,spatgev_data$X,include_tau=include_tau)
  }
  
  return(list(linpred=linpred,X=spatgev_data$X,X_unscaled=X_unscaled,X_clim=X_clim,X_clim_unscaled=X_clim_unscaled))
}



calculate_parameters_clim=function(linpred_res,targlon,targlat){
  X_clim=linpred_res$X_clim_unscaled
  uniquelocs=unique(data.table(X_clim)[,.(lon,lat)])
  ind=which(X_clim[,lon]==targlon & X_clim[,lat]==targlat)
  
  my_linpred=linpred_res$linpred
  
  mu=unlist(my_linpred$mu[ind])
  kappa=unlist(my_linpred$kappa[ind])
  xi=unlist(my_linpred$xi[ind])
  
  return(list(mu=mu,kappa=kappa,xi=xi,lon=targlon,lat=targlat))
}

#---------------------------------------------------------------------------------------------------#
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")

load("/nr/project/stat/Impetus4Change/Res/cv_temporal_subset_spacetau/nonspat1_allyears.Rdata")
mcmc_temp=mcmc_all
linpred_temp=summarize_results(mcmc_temp,amax_data,fit_type="annual",include_tau=1)
linpred_temp0=summarize_results(mcmc_temp,amax_data,fit_type="annual",clim_res=TRUE,include_tau=1)

#---------------------------------------------------------------------------------------------------#
load("/nr/project/stat/Impetus4Change/Res/cv_spatial/nonspat1_allyears.Rdata")
mcmc_spat=mcmc_all
linpred_spat=summarize_results(mcmc_spat,amax_data,fit_type="full_spatial",include_tau=1)
#---------------------------------------------------------------------------------------------------#

for(locind in 1:12){
  uniquelocs=unique(linpred_temp0$X_clim_unscaled[,.(lon,lat)])
  targlon=uniquelocs$lon[locind]
  targlat=uniquelocs$lat[locind]
  
  #---------Histograms of mu, kappa, xi for this location------------
  #temporal:
  params_temp=calculate_parameters_clim(linpred_temp0,targlon=targlon,targlat=targlat)
  mu_temp=params_temp$mu
  kappa_temp=params_temp$kappa
  xi_temp=params_temp$xi
  
  #spatial:
  params_spat=calculate_parameters_clim(linpred_spat,targlon=targlon,targlat=targlat)
  mu_spat=params_spat$mu
  kappa_spat=params_spat$kappa
  xi_spat=params_spat$xi
  
  par(mfrow=c(2,3))
  truehist(unlist(mu_temp),xlim=c(0,25),ylim=c(0,2))
  abline(v=median(unlist(mu_temp)),col="blue")
  abline(v=median(unlist(mu_spat)),col="red")
  abline(v=freqres[locind,1],col="orange")
  title(locind)
  
  
  truehist(unlist(kappa_temp),xlim=c(0,10),ylim=c(0,10))
  abline(v=median(unlist(kappa_temp)),col="blue")
  abline(v=median(unlist(kappa_spat)),col="red")
  abline(v=freqres[locind,2],col="orange")
  title(locind)
  
  truehist(unlist(xi_temp),xlim=c(-3,3),ylim=c(0,10))
  abline(v=median(unlist(xi_temp)),col="blue")
  abline(v=median(unlist(xi_spat)),col="red")
  abline(v=freqres[locind,3],col="orange")
  title(locind)
  
  truehist(unlist(mu_spat),xlim=c(0,25),ylim=c(0,2))
  abline(v=median(unlist(mu_temp)),col="blue")
  abline(v=median(unlist(mu_spat)),col="red")
  abline(v=freqres[locind,1],col="orange")
  title(locind)
  
  truehist(unlist(kappa_spat),xlim=c(0,10),ylim=c(0,10))
  abline(v=median(unlist(kappa_temp)),col="blue")
  abline(v=median(unlist(kappa_spat)),col="red")
  abline(v=freqres[locind,2],col="orange")
  title(locind)
  
  truehist(unlist(xi_spat),xlim=c(-3,3),ylim=c(0,10))
  abline(v=median(unlist(xi_temp)),col="blue")
  abline(v=median(unlist(xi_spat)),col="red")
  abline(v=freqres[locind,3],col="orange")
  title(locind)
  
}

#--------------------------------------------------------------#
# Look at GEV-curves for individual locations#

#locind=1
par(mfrow=c(1,3),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
for(locind in 1:12){
  targlon=uniquelocs$lon[locind]; targlat=uniquelocs$lat[locind]
  
  #temporal:
  params_temp=calculate_parameters_clim(linpred_temp0,targlon=targlon,targlat=targlat)
  params_temp=cbind(params_temp$mu,params_temp$kappa,params_temp$xi)
  params_temp[which(params_temp[,2]<0),2]=0.01

  #spatial:
  params_spat=calculate_parameters_clim(linpred_spat,targlon=targlon,targlat=targlat)
  params_spat=cbind(params_spat$mu,params_spat$kappa,params_spat$xi)
  params_spat[which(params_spat[,2]<0),2]=0.01
  
  
  
  #---Gev distributions----------#
  this_retper=c(2,5,10,20,25,30,40,50,60,70,80,90,100)
  retper_temp=t(apply(params_spat[1:90001,],1,gevwrap,retper=this_retper))
  retper_spat=t(apply(params_temp[1:90001,],1,gevwrap,retper=this_retper))
  retper_freq=gevwrap(freqres[locind,],retper=this_retper)
  
  retper_median_temp=apply(retper_temp,2,quantile,probs=c(0.05,0.5,0.95))
  retper_median_spat=apply(retper_spat,2,quantile,probs=c(0.05,0.5,0.95))
  
  plot(this_retper,retper_median_temp[2,],type="o",col="darkblue",ylim=c(0,100),xlab="Return period",ylab="Return level (mm)");grid()
  lines(this_retper,retper_median_temp[1,],lty=2,col="darkblue")
  lines(this_retper,retper_median_temp[3,],lty=2,col="darkblue")
  lines(this_retper,retper_freq,col="orange",type="o")
  
  lines(this_retper,retper_median_spat[2,],lty=1,col="red",type="o")
  lines(this_retper,retper_median_spat[1,],lty=2,col="red")
  lines(this_retper,retper_median_spat[3,],lty=2,col="red")
  title(locind)

}

#------------------------------#
prob_inc_mu=round(rbind(colMeans(mcmc_temp$M[,,1]),colMeans(mcmc_spat$M[,,1])),2)
prob_inc_kappa=round(rbind(colMeans(mcmc_temp$M[,,2]),colMeans(mcmc_spat$M[,,2])),2)
prob_inc_xi=round(rbind(colMeans(mcmc_temp$M[,,3]),colMeans(mcmc_spat$M[,,3])),2)

coef_mu=rbind(apply(mcmc_temp$THETA[,,1],2,median),apply(mcmc_spat$THETA[,,1],2,median))
coef_kappa=rbind(apply(mcmc_temp$THETA[,,2],2,median),apply(mcmc_spat$THETA[,,2],2,median))
coef_xi=rbind(apply(mcmc_temp$THETA[,,3],2,median),apply(mcmc_spat$THETA[,,3],2,median))

colnames(prob_inc_mu)=colnames(prob_inc_kappa)=colnames(prob_inc_xi)=mcmc_temp$covariates
colnames(coef_mu)=colnames(coef_kappa)=colnames(coef_xi)=mcmc_temp$covariates


coef_mu
coef_kappa
coef_xi
#------------------------------#
