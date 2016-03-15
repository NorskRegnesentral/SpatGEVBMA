##### JUST INSTALLATION STUFF HERE ####

# Specific for MJ

### Install new STAR version of SpatialGEVBMA by running this command:
# cd /nr/samba/user/jullum/Prosjekter/FlomQ/Git_SpatialGEVBMA
# R CMD INSTALL --library="/nr/samba/user/jullum/Rlibstest/new_packages" SpatialGEVBMA

### Install OLD version of SpatialGEVBMA by running this command:
# cd /nr/samba/user/jullum/Prosjekter/FlomQ/SpatialGEVBMA_package_version_160202
# R CMD INSTALL --library="/nr/samba/user/jullum/Rlibstest/old_packages" SpatialGEVBMA

### Install STAR version of spatial.gev.bma by running this command:
# cd /nr/samba/user/jullum/Prosjekter/FlomQ/Old_spatial.gev.bma_star
# R CMD INSTALL --library="/nr/samba/user/jullum/Rlibstest/test_packages" spatial.gev.bma


#############################################



rm(list = ls())

### New version of SpatialGEVBMA
library("SpatialGEVBMA", lib.loc="/nr/samba/user/jullum/Rlibstest/new_packages")

data(norway)

S=norway$S
X=norway$X
Y.list=norway$Y.list


set.seed(123)

Rprof()
a <- spatial.gev.bma(Y.list, X, S, 10^5,print.every=2000,nonspatial = TRUE,log.kappa=TRUE)
Rprof(NULL)
a.summary=summaryRprof()

a.summary$sampling.time

detach("package:SpatialGEVBMA",unload=TRUE)



### For comparison with spatial.gev.bma 
library(spatial.gev.bma,lib.loc="/nr/samba/user/jullum/Rlibstest/test_packages")

set.seed(123)

Rprof()
b <- spatial.gev.bma(Y.list, X, S, 10^5,print.every=2000,nonspatial = TRUE,log.kappa=TRUE)
Rprof(NULL)
b.summary=summaryRprof()

b.summary$sampling.time
detach("package:spatial.gev.bma",unload=TRUE)





##############


colsd=function(x)
{
  apply(x, 2, sd)
}

means.THETA.a=NULL
means.M.a=NULL
means.TAU.a=NULL
means.ACCEPT.TAU.a=NULL
sd.THETA.a=NULL
sd.M.a=NULL
sd.TAU.a=NULL
sd.ACCEPT.TAU.a=NULL


for (j in 1:3)
{
  
  means.THETA.a=c(means.THETA.a,colMeans(a$THETA[,,j]))
  means.M.a=c(means.M.a,colMeans(a$M[,,j]))
  means.TAU.a=c(means.TAU.a,colMeans(a$TAU[,,j]))
  means.ACCEPT.TAU.a=c(means.ACCEPT.TAU.a,colMeans(a$ACCEPT.TAU[,,j]))
  
  sd.THETA.a=c(sd.THETA.a,colsd(a$THETA[,,j]))
  sd.M.a=c(sd.M.a,colsd(a$M[,,j]))
  sd.TAU.a=c(sd.TAU.a,colsd(a$TAU[,,j]))
  sd.ACCEPT.TAU.a=c(sd.ACCEPT.TAU.a,colsd(a$ACCEPT.TAU[,,j]))
  
}
means.LAMBDA.a=colMeans(a$LAMBDA)
means.ALPHA.a=colMeans(a$ALPHA)

sd.LAMBDA.a=colsd(a$LAMBDA)
sd.ALPHA.a=colsd(a$ALPHA)


means.THETA.b=NULL
means.M.b=NULL
means.TAU.b=NULL
means.ACCEPT.TAU.b=NULL

sd.THETA.b=NULL
sd.M.b=NULL
sd.TAU.b=NULL
sd.ACCEPT.TAU.b=NULL


for (j in 1:3)
{
  means.THETA.b=c(means.THETA.b,colMeans(b$THETA[,,j]))
  means.M.b=c(means.M.b,colMeans(b$M[,,j]))
  means.TAU.b=c(means.TAU.b,colMeans(b$TAU[,,j]))
  means.ACCEPT.TAU.b=c(means.ACCEPT.TAU.b,colMeans(b$ACCEPT.TAU[,,j]))
  
  sd.THETA.b=c(sd.THETA.b,colsd(b$THETA[,,j]))
  sd.M.b=c(sd.M.b,colsd(b$M[,,j]))
  sd.TAU.b=c(sd.TAU.b,colsd(b$TAU[,,j]))
  sd.ACCEPT.TAU.b=c(sd.ACCEPT.TAU.b,colsd(b$ACCEPT.TAU[,,j]))
}
means.LAMBDA.b=colMeans(b$LAMBDA)
means.ALPHA.b=colMeans(b$ALPHA)

sd.LAMBDA.b=colsd(b$LAMBDA)
sd.ALPHA.b=colsd(b$ALPHA)

means.THETA.a=means.THETA.a
means.THETA.b=means.THETA.b

#setwd("/nr/samba/user/jullum/Prosjekter/FlomQ/")
pdf(file="comp_with_spatial.gev.bma.star3.pdf")
plot(abs(means.THETA.a),log='y',ylim=range(abs(c(means.THETA.a,means.THETA.b,means.THETA.a+sd.THETA.a,means.THETA.b+sd.THETA.b))))
points(1:length(means.THETA.a)+0.2,abs(means.THETA.b),col=2)
segments(1:length(means.THETA.a),abs(means.THETA.a),1:length(means.THETA.a),abs(means.THETA.a)+sd.THETA.a)
segments((1:length(means.THETA.b))+0.2,abs(means.THETA.b),1:length(means.THETA.a)+0.2,abs(means.THETA.b)+sd.THETA.b,col=2)

plot(abs(means.M.a),log='y',ylim=range(abs(c(means.M.a,means.M.b,means.M.a+sd.M.a,means.M.b+sd.M.b))))
points(1:length(means.M.a)+0.2,abs(means.M.b),col=2)
segments(1:length(means.M.a),abs(means.M.a),1:length(means.M.a),abs(means.M.a)+sd.M.a)
segments((1:length(means.M.b))+0.2,abs(means.M.b),1:length(means.M.a)+0.2,abs(means.M.b)+sd.M.b,col=2)

par(cex=0.5)
plot(abs(means.TAU.a),log='y',ylim=range(abs(c(means.TAU.a,means.TAU.b,means.TAU.a+sd.TAU.a,means.TAU.b+sd.TAU.b))))
points(1:length(means.TAU.a)+0.2,abs(means.TAU.b),col=2)
segments(1:length(means.TAU.a),abs(means.TAU.a),1:length(means.TAU.a),abs(means.TAU.a)+sd.TAU.a)
segments((1:length(means.TAU.b))+0.2,abs(means.TAU.b),1:length(means.TAU.a)+0.2,abs(means.TAU.b)+sd.TAU.b,col=2)

plot(abs(means.ACCEPT.TAU.a),log='y',ylim=range(abs(c(means.ACCEPT.TAU.a,means.ACCEPT.TAU.b,means.ACCEPT.TAU.a+sd.ACCEPT.TAU.a,means.ACCEPT.TAU.b+sd.ACCEPT.TAU.b))))
points(1:length(means.ACCEPT.TAU.a)+0.2,abs(means.ACCEPT.TAU.b),col=2)
segments(1:length(means.ACCEPT.TAU.a),abs(means.ACCEPT.TAU.a),1:length(means.ACCEPT.TAU.a),abs(means.ACCEPT.TAU.a)+sd.ACCEPT.TAU.a)
segments((1:length(means.ACCEPT.TAU.b))+0.2,abs(means.ACCEPT.TAU.b),1:length(means.ACCEPT.TAU.a)+0.2,abs(means.ACCEPT.TAU.b)+sd.ACCEPT.TAU.b,col=2)

par(cex=1)
plot(abs(means.LAMBDA.a),log='y',ylim=range(abs(c(means.LAMBDA.a,means.LAMBDA.b,means.LAMBDA.a+sd.LAMBDA.a,means.LAMBDA.b+sd.LAMBDA.b))))
points(1:length(means.LAMBDA.a)+0.02,abs(means.LAMBDA.b),col=2)
segments(1:length(means.LAMBDA.a),abs(means.LAMBDA.a),1:length(means.LAMBDA.a),abs(means.LAMBDA.a)+sd.LAMBDA.a)
segments((1:length(means.LAMBDA.b))+0.02,abs(means.LAMBDA.b),1:length(means.LAMBDA.a)+0.02,abs(means.LAMBDA.b)+sd.LAMBDA.b,col=2)

plot(abs(means.ALPHA.a),log='y',ylim=range(abs(c(means.ALPHA.a,means.ALPHA.b,means.ALPHA.a+sd.ALPHA.a,means.ALPHA.b+sd.ALPHA.b))))
points(1:length(means.ALPHA.a)+0.02,abs(means.ALPHA.b),col=2)
segments(1:length(means.ALPHA.a),abs(means.ALPHA.a),1:length(means.ALPHA.a),abs(means.ALPHA.a)+sd.ALPHA.a)
segments((1:length(means.ALPHA.b))+0.02,abs(means.ALPHA.b),1:length(means.ALPHA.a)+0.02,abs(means.ALPHA.b)+sd.ALPHA.b,col=2)

dev.off()


