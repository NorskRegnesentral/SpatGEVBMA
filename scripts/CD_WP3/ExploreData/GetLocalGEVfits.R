rm(list=ls())

library(readxl) ## To read in data set 
library(fExtremes) ## To fit local GEV distributions 
library(data.table)



datapath = "~/Documents/Data/ExtremePrecipNorway/"
my.datafile = paste(datapath, "AM_all_2023.xlsx", sep="")

durations = c(10,15,20,30,45,60,90,120,180,360,720,1440)
DT <- data.table("station" = character(),    
                    "duration" = double(),
                    "mu" = double(),
                    "sigma" = double(),
                    "xi" = double())

for(d in durations) {
  data = read_excel(my.datafile, sheet=paste(d,"min",sep=""))         
  N = dim(data)[2] 
  for(i in 2:N) {
    y = as.numeric(data[[i]])
    ind = which(is.na(y))
    y = y[-ind]
    pwm.coefs <- gevFit(y, type="pwm")@fit$par.ests
    tmp = list(station=names(data)[i],
               duration=d,
               mu=pwm.coefs[2],
               sigma=pwm.coefs[3],
               xi=pwm.coefs[1])
    DT = rbind(DT, tmp)
  }
}

save(DT,file="LocalGEVfits.RData")
