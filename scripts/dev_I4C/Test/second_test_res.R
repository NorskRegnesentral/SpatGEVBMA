library(data.table)
source("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/R/gev.R")
source("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/R/temporal_spatgev.R")
source("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/R/RcppExports.R")

makeLinPred=function(R,burnin=1e2,X,nonspatial=FALSE){
  #kappa could be logged.
  
  n=dim(R$TAU)[1]
  mu=list(); kappa=list(); xi=list()
  
  #if(nonspatial==FALSE){
    for(stedind in 1:dim(X)[1]){
      #Param nr 1 (mu?):
      mu[[stedind]]=rowSums((R$THETA[burnin:n,,1])*matrix(rep(X[stedind,],length(burnin:n)),byrow=TRUE,ncol=length(X[stedind,])))+
        R$TAU[burnin:n,stedind,1]
      
      kappa[[stedind]]=rowSums((R$THETA[burnin:n,,2])*matrix(rep(X[stedind,],length(burnin:n)),byrow=TRUE,ncol=length(X[stedind,])))+
        R$TAU[burnin:n,stedind,2]
      
      
      xi[[stedind]]=rowSums((R$THETA[burnin:n,,3])*matrix(rep(X[stedind,],length(burnin:n)),byrow=TRUE,ncol=length(X[stedind,])))+
        R$TAU[burnin:n,stedind,3]
    }
    
    if(R$log.kappa==TRUE){
      kappa=lapply(kappa,exp)
    }
 # }else{
    
    
    
 # }
  
  
  
  
  return(list(kappa=kappa,mu=mu,xi=xi))
}

#--------------------------------------------------------------------------------------------------------------------------------#
amax_data=fread(file="/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/dev_I4C/Data/AM60_cov.csv")[,.(lon,lat,year,masl,stid,
                                                                                                       wetterdays_MAM,wetterdays_JJA,wetterdays_SON,wetterdays_annual,
                                                                                                       precip_MAM,precip_JJA,precip_SON,precip_annual,
                                                                                                       temp_MAM,temp_JJA,temp_SON,temp_annual,y)]

amax_data=amax_data[lon<7&lat<62 & lat>58]

spatgev_data=make_temporal_spatgev_data(amax_data,TRUE)

load("/nr/project/stat/Impetus4Change/Res/mcmc_temporal_spatial50000_hourly.Rdata")

burnin=5000
n=dim(mcmc_res$THETA)[1]

#------------------------------------------------------------------------------------------------------------------------"
linpred=makeLinPred(mcmc_res,spatgev_data$X,burnin=burnin)

Modelmat=rbind(colMeans(mcmc_res$M[burnin:n,,1]),colMeans(mcmc_res$M[burnin:n,,2]),colMeans(mcmc_res$M[burnin:n,,3]))*100
colnames(Modelmat)=colnames(spatgev_data$X)

coefs=round(rbind(apply(mcmc_res$THETA[,,1],2,median),
            apply(mcmc_res$THETA[,,2],2,median),
            apply(mcmc_res$THETA[,,3],2,median)),3)

colnames(coefs)=colnames(spatgev_data$X)

print(round(Modelmat,1))
print(coefs)
#------------------------------------------------------------------------------------------------------------------------"


yrloc=amax_data[,.(stid,year)]
yrloc[,num:=1:(dim(amax_data)[1])]
uniquelocs=unique(yrloc$stid)

#if non-spatial=TRUE:
plot(NA,ylim=c(-100,100),xlim=c(1960,2022))
for(j in 1:length(uniquelocs)){
  currloc=yrloc[stid==uniquelocs[j]]
  curr_postmedian=apply(mcmc_res$TAU[,currloc$num,3],2,mean)
  
  lines(currloc$year,curr_postmedian,type="o",col=j,lty=j);grid()

}


#if non-spatial = FALSE:
for(j in 1:length(uniquelocs)){
  currloc=yrloc[stid==stid_numyr]
  curr_postmedian=apply(mcmc_res$TAU[,currloc$num,1],2,mean)
  
  lines(currloc$year,curr_postmedian,type="o",col=j,lty=j);grid()
  
}


#Alle får samme tidseffekt --> Bra.

colMeans(mcmc_res$ACCEPT.TAU[,,1])
