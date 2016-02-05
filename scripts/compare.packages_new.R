### Install old version of package by running this command:

# cd /nr/samba/user/jullum/Prosjekter/FlomQ/SpatialGEVBMA_package_version_160202
# R CMD INSTALL --library="/home/martin/R/old_packages" SpatialGEVBMA


### Start by installing the package

## Run R CMD INSTALL SpatialGEVBMA from terminal when being in the 
## folder above the one where the package is located

## For Martin that is:

# cd /nr/samba/user/jullum/Prosjekter/FlomQ/Git_SpatialGEVBMA/
# R CMD INSTALL SpatialGEVBMA


rm(list = ls())

library(SpatialGEVBMA)
#install.packages("spatial.gev.bma")
#library("spatial.gev.bma")

data(norway)
#attach(norway)

#S=norway$S
#X=norway$X
#Y.list=norway$Y.list

S=norway$S
X=norway$X
Y.list=norway$Y.list


set.seed(123)

Rprof()
a <- SpatialGEVBMA::spatial.gev.bma(Y.list, X, S, 10^3,print.every=200)
Rprof(NULL)
a.summary=summaryRprof()

a.summary$sampling.time

### For comparison with the old version of the package 
library(SpatialGEVBMA,lib.loc="/home/martin/R/old_packages")


set.seed(123)

Rprof()
b <- spatial.gev.bma::spatial.gev.bma(Y.list, X, S, 10^3,print.every=200)
Rprof(NULL)
b.summary=summaryRprof()

b.summary$sampling.time



###


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


pdf(file="test.pdf")
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





##### OLD stuff underneath








plot(density(a$THETA[,8,3]))
lines(density(b$THETA[,8,3]),col=2)


sum(a$M==0)/length(a$M)
sum(b$M==0)/length(a$M)  
  
plot(density(a$TAU))
lines(density(b$TAU),col=2,lty=2)

plot(density(a$LAMBDA))
lines(density(b$LAMBDA),col=2,lty=2)

plot(density(a$ALPHA))
lines(density(b$ALPHA),col=2,lty=2)

sum(a$ACCEPT.TAU==0)/length(a$M)
sum(b$ACCEPT.TAU==0)/length(a$M)  

plot(density(a$S))
lines(density(b$S),col=2,lty=2)


a$TAU[190:195,8,2]
b$TAU[190:195,8,2]


#> a$TAU[131,1,3]
#[1] -0.0009564286
#> b$TAU[131,1,3]
#[1] -0.0006841234
#> a$TAU[131,1,3]/b$TAU[131,1,3]
#[1] 1.398035



# Forskjelleige
a$THETA[131,1,3]
b$THETA[131,1,3]

# Like
a$M[131,1,3]
b$M[131,1,3]

# Forskjellige
a$TAU[131,1,3]
b$TAU[131,1,3]

# Forskjellige
a$LAMBDA[131,3]
b$LAMBDA[131,3]

# Forskjellige
a$ALPHA[131,3]
b$ALPHA[131,3]

# Like
a$ACCEPT.TAU[131,1,3]
b$ACCEPT.TAU[131,1,3]




#> a.summary$sampling.time
#[1] 20.16
#> b.summary$sampling.time
#[1] 19.9

#> a.summary$sampling.time
#[1] 18.46
#> b.summary$sampling.time
#[1] 20.16
