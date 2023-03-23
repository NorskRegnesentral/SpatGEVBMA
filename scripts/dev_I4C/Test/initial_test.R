library(SpatGEVBMA);
library(data.table)
setwd("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA")
source("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/R/gev.R")

#Sys.setenv(OPENBLAS_NUM_THREADS=1)



makeLinPred=function(R,burnin=1e2,X){
  #kappa could be logged.
  
  n=dim(R$TAU)[1]
  
  mu=list(); kappa=list(); xi=list()
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
  
  return(list(kappa=kappa,mu=mu,xi=xi))
}


#----Upload data-----#

load(paste("data/norway.RData"))

Y=norway$Y.list
X=norway$X
S=norway$S

X.list=X;Y.list=Y;S=S


#----Pretend that we only have one location-----#
Yoneloc=list()
for(j in 1:length(Y[[41]])){
  Yoneloc[[j]]=(Y[[41]])[j]
}

targlon=data.table(X)[41,lon];targlat=data.table(X)[41,lat]
Xoneloc=copy(data.table(X[1:length(Y[[41]]),]))
Xoneloc[,lon:=targlon];Xoneloc[,lat:=targlat]
Xoneloc=Xoneloc[,.(lon,lat)];Xoneloc[,year:=1950:(1950+length(Y[[41]])-1)]
Xoneloc=cbind(intercept=1,Xoneloc[,.(year)])
Soneloc=cbind(0,Xoneloc$year) #representerer år.

Xoneloc=as.matrix(data.frame(Xoneloc))
Soneloc=as.matrix(data.frame(Soneloc))

#--------Prior settings------------#
p <- dim(X)[2]
prior <- NULL
prior$mu$alpha.a <- 2
prior$mu$alpha.b <- 6
prior$mu$lambda.a <- 2
prior$mu$lambda.b <- 2

prior$kappa$alpha.a <- 2
prior$kappa$alpha.b <- 2
prior$kappa$lambda.a <- 1.5
prior$kappa$lambda.b <- 1.5

prior$xi$alpha.a <- 2
prior$xi$alpha.b <- 1
prior$xi$lambda.a <- 2
prior$xi$lambda.b <- 1

prior$mu$beta.0 <- c(8,rep(0, p-1))
#---------------------------------#

n.reps <- 5000#100 #2e5
#if(0){
  Roneloc <- spatial.gev.bma(Yoneloc, Xoneloc, Soneloc,
                             n.reps, prior, print.every = 100,
                             i4c=TRUE)
  
  #save(Roneloc,file="Test2023/res/oneloc_full.RData")
  
  #R <- spatial.gev.bma(Y, X, S, n.reps, prior, print.every = 1000,i4c=FALSE)
  #save(R,file="Test2023/res/orig_full.RData")
  
 
#}
  
  if(0){
  
  load("Test2023/res/oneloc_full.RData")
  linpred_oneloc=makeLinPred(Roneloc,burnin=20000,Xoneloc)
  
  #load("Test2023/res/oneloc100.RData")
  #linpred_oneloc=makeLinPred(Roneloc,burnin=20,Xoneloc)
  
  burnin=20000
  n=dim(Roneloc$THETA)[1]
  mu_intercept=median(Roneloc$THETA[burnin:n,1,1])
  kappa_intercept=median(Roneloc$THETA[burnin:n,1,2])
  xi_intercept=median(Roneloc$THETA[burnin:n,1,3])
  
  Ydata=unlist(Yoneloc)
  par(mfrow=c(3,3))
  for(yearnum in 1:9){
    hist(linpred_oneloc$mu[[yearnum]]);abline(v=mu_intercept,col="red");abline(v=Ydata[yearnum],col="green")
    hist(linpred_oneloc$kappa[[yearnum]]);abline(v=kappa_intercept,col="red")
    hist(linpred_oneloc$xi[[yearnum]]);abline(v=xi_intercept,col="red")
  }
  
  inclusion_prob=round(rbind(colMeans(Roneloc$M[burnin:n,,1]),
        colMeans(Roneloc$M[burnin:n,,2]),
        colMeans(Roneloc$M[burnin:n,,3]))*100,2)
  
  colnames(inclusion_prob)=c("Intercept","Year cov")
  rownames(inclusion_prob)=c("mu","kappa","xi")
  inclusion_prob
  
  #Årseffekt:
  dim(Roneloc$TAU[,,1])
  dim(Roneloc$TAU[burnin:n,,1])
  
  time_effect_mu=apply(Roneloc$TAU[burnin:n,,1],2,quantile,probs=c(0.05,0.5,0.95))
  time_effect_kappa=apply(Roneloc$TAU[burnin:n,,2],2,quantile,probs=c(0.05,0.5,0.95))
  time_effect_xi=apply(Roneloc$TAU[burnin:n,,3],2,quantile,probs=c(0.05,0.5,0.95))
  
  par(mfrow=c(3,1))
  plot(Xoneloc[,2],time_effect_mu[2,],type="p",xlab="Year",ylim=c(-20,20));grid()
  lines(Xoneloc[,2],time_effect_mu[1,],lty=2)
  lines(Xoneloc[,2],time_effect_mu[3,],lty=2)
  lines(Xoneloc[,2],unlist(Yoneloc),lty=1,col="green")
  title("mu")
  
  plot(Xoneloc[,2],time_effect_kappa[2,],type="p",xlab="Year",ylim=c(-2,2));grid()
  lines(Xoneloc[,2],time_effect_kappa[1,],lty=2)
  lines(Xoneloc[,2],time_effect_kappa[3,],lty=2)
  title("kappa")
  
  plot(Xoneloc[,2],time_effect_xi[2,],type="p",xlab="Year",ylim=c(-2,2));grid()
  lines(Xoneloc[,2],time_effect_xi[1,],lty=2)
  lines(Xoneloc[,2],time_effect_xi[3,],lty=2)
  title("xi")
  
  #--------------Test with two locations-----------------------------#
  #(We don't know which years we have data):
  Ytwoloc=list();k=1
  for(j in 1:length(Y[[41]])){
    Ytwoloc[[k]]=(Y[[41]])[j]
    k=k+1
  }
  
  
  for(j in 1:length(Y[[40]])){
    Ytwoloc[[k]]=(Y[[40]])[j]
    k=k+1
  }
  
  
  
  targlon=data.table(X)[41,lon];targlat=data.table(X)[41,lat]
  Xoneloc=copy(data.table(X[1:length(Y[[41]]),]))
  Xoneloc[,lon:=targlon];Xoneloc[,lat:=targlat]
  Xoneloc=Xoneloc[,.(lon,lat)];Xoneloc[,year:=1950:(1950+length(Y[[41]])-1)]
  Xoneloc=cbind(intercept=1,Xoneloc[,.(year)],lon=targlon,lat=targlat)
  
  firstloc=Xoneloc
  
  targlon=data.table(X)[40,lon];targlat=data.table(X)[40,lat]
  Xoneloc=copy(data.table(X[1:length(Y[[40]]),]))
  Xoneloc[,lon:=targlon];Xoneloc[,lat:=targlat]
  Xoneloc=Xoneloc[,.(lon,lat)];Xoneloc[,year:=tail(firstloc$year,length(Y[[40]]))]
  Xoneloc=cbind(intercept=1,Xoneloc[,.(year)],lon=targlon,lat=targlat)
  
  Xtwoloc=rbind(firstloc,Xoneloc)
  
  Stwoloc=cbind(0,Xtwoloc$year) #representerer år.
  
  Xtwoloc=as.matrix(data.frame(Xtwoloc))
  Stwoloc=as.matrix(data.frame(Stwoloc))
  
  n.reps <- 2000#100 #2e5
  if(0){
    Rtwoloc <- spatial.gev.bma(Ytwoloc, Xtwoloc, Stwoloc,
                               n.reps, prior, print.every = 100,
                               i4c=TRUE)
    
    #save(Rtwoloc,file="Test2023/res/twoloc_full.RData")
  
  }
  
  
  load("Test2023/res/twoloc_full.RData")
  Roneloc=Rtwoloc
  linpred_oneloc=makeLinPred(Roneloc,burnin=100,Xtwoloc)
  }
  
  
  #----------------------------
  #go to folder
  #----------------------------
  # export OPENBLAS_NUM_THREADS=1
  #R CMD BATCH initial_tets.R
  #----------------------------
  #bach.c file.
  #https://github.com/SeasonalForecastingEngine/NordpoolDemand/blob/master/R/pca_analysis.R