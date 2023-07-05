library(ggplot2)
return_periods=c(2,5,10,25,50,100,200)
quantiles=1-1/return_periods

#Upload hourly data (change 60 to 1440 to get daily data):
amax_data=fread(file="scripts/dev_I4C/Data/AM60_cov.csv")[,.(year,masl,stid,lon,lat,
                                                             wetterdays_JJA,wetterdays_annual,
                                                             precip_JJA,precip_annual,
                                                             temp_JJA,temp_annual,year,y)]


res=fread("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/bGev_tutorial/inlares.csv")
cor(res[year_specific==TRUE & spatial==TRUE,mu],res[year_specific==TRUE & spatial==TRUE,y])
cor(res[year_specific==FALSE & spatial==TRUE,mu],res[year_specific==FALSE & spatial==TRUE,y])

uniquelocs=unique(res[,.(lon,lat)])

p=1
allres=list()
for(locnum in 1:(dim(uniquelocs)[1])){
  print(locnum)
  thisloc=res[lon==uniquelocs$lon[locnum] & lat==uniquelocs$lat[locnum]]
  years=unique(thisloc[,year])
  
  for(yy in years){
    spatgevres_baseline=c()
    spatgevres_yearspecific=c()
    gevres=c()
    thislocyear=thisloc[year==yy & spatial==TRUE]
    for(j in 1:length(quantiles)){
      spatgevres_baseline[j]=qbgev(quantiles[j],mu=thislocyear[year_specific==FALSE,mu],
                       sigma=thislocyear[year_specific==FALSE,sigma],
                       xi=thislocyear[year_specific==FALSE,xi],p_b=0.2,s=5)
      
      spatgevres_yearspecific[j]=qbgev(quantiles[j],mu=thislocyear[year_specific==TRUE,mu],
                                       sigma=thislocyear[year_specific==TRUE,sigma],
                                       xi=thislocyear[year_specific==TRUE,xi],p_b=0.2,s=5)
      
      gevres[j]=evd::qgev(quantiles[j],thislocyear[year_specific==FALSE,loc_gev],
                          thislocyear[year_specific==FALSE,scale_gev],thislocyear[year_specific==FALSE,shape_gev])
    }
    
    spatgevres_yearspecific=spatgevres_yearspecific*thislocyear[year_specific==TRUE,clim_sd]+thislocyear[year_specific==TRUE,clim_mean]
    spatgevres_baseline=spatgevres_baseline*thislocyear[year_specific==FALSE,clim_sd]+thislocyear[year_specific==FALSE,clim_mean]
  
    spatgevres_yearspecific=data.table(return_level=spatgevres_yearspecific,return_periods,year=yy,y=thislocyear[1,y],lon=uniquelocs$lon[locnum],lat=uniquelocs$lat[locnum],type="bgev_year")
    spatgevres_baseline=data.table(return_level=spatgevres_baseline,return_periods,year=yy,y=thislocyear[1,y],lon=uniquelocs$lon[locnum],lat=uniquelocs$lat[locnum],type="bgev_base")
    gevres=data.table(return_level=gevres,return_periods,year=yy,y=thislocyear[1,y],lon=uniquelocs$lon[locnum],lat=uniquelocs$lat[locnum],type="gev")
   
    allres[[p]]=rbind(spatgevres_yearspecific,spatgevres_baseline,gevres)
    p=p+1
  }
}

allres=rbindlist(allres)
#fwrite(allres,"/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/bGev_tutorial/retlevel_res.csv")

allres=fread("/nr/samba/user/roksvag/GitRepo/SpatGEVBMA/scripts/bGev_tutorial/retlevel_res.csv")

library(SpatGEVBMA)
par(mfrow=c(1,3))
plot(allres[return_periods==2 & type=="bgev_year",y,return_level],xlim=c(0,20));lines(c(0,1000),c(0,1000),col="red")
plot(allres[return_periods==2 & type=="bgev_base",y,return_level],xlim=c(0,20));lines(c(0,1000),c(0,1000),col="red")
plot(allres[return_periods==2 & type=="gev",y,return_level],xlim=c(0,20));lines(c(0,1000),c(0,1000),col="red")

par(mfrow=c(1,4))
plot(allres[return_periods==2 & type=="bgev_base",return_level],
     allres[return_periods==2 & type=="bgev_year",return_level]);lines(c(0,1000),c(0,1000),col="red")

plot(allres[return_periods==10 & type=="bgev_base",return_level],
     allres[return_periods==10 & type=="bgev_year",return_level]);lines(c(0,1000),c(0,1000),col="red")

plot(allres[return_periods==50 & type=="bgev_base",return_level],
     allres[return_periods==50 & type=="bgev_year",return_level]);lines(c(0,1000),c(0,1000),col="red")

plot(allres[return_periods==100 & type=="bgev_base",return_level],
     allres[return_periods==100 & type=="bgev_year",return_level]);lines(c(0,1000),c(0,1000),col="red")



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

coefficients=unique(res[spatial==TRUE ,.(year,year_specific,beta0.5_intercept,beta0.5_lat,beta0.5_lon,beta0.5_masl,beta0.5_precip_JJA,beta0.5_temp_JJA,beta0.5_wetterdays_JJA,
                                         beta0.975_intercept,beta0.975_lat,beta0.975_lon,beta0.975_masl,beta0.975_precip_JJA,beta0.975_temp_JJA,beta0.975_wetterdays_JJA,
                                         beta0.025_intercept,beta0.025_lat,beta0.025_lon,beta0.025_masl,beta0.025_precip_JJA,beta0.025_temp_JJA,beta0.025_wetterdays_JJA)])


pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//precip_JJA.pdf",width = 8, height = 4)
ggplot(coefficients,aes(x=year,y=beta0.5_precip_JJA,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=beta0.025_precip_JJA, ymax=beta0.975_precip_JJA))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//temp_JJA.pdf",width = 8, height = 4)
ggplot(coefficients,aes(x=year,y=beta0.5_temp_JJA,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=beta0.025_temp_JJA, ymax=beta0.975_temp_JJA))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//wetterdays_JJA.pdf",width = 8, height = 4)
ggplot(coefficients,aes(x=year,y=beta0.5_wetterdays_JJA,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=beta0.025_wetterdays_JJA, ymax=beta0.975_wetterdays_JJA))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//lon.pdf",width = 8, height = 4)
ggplot(coefficients,aes(x=year,y=beta0.5_lon,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=beta0.025_lon, ymax=beta0.975_lon))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//lat.pdf",width = 8, height = 4)
ggplot(coefficients,aes(x=year,y=beta0.5_lat,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=beta0.025_lat, ymax=beta0.975_lat))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//masl.pdf",width = 8, height = 4)
ggplot(coefficients,aes(x=year,y=beta0.5_masl,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=beta0.025_masl, ymax=beta0.975_masl))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//intercept.pdf",width = 8, height = 4)
ggplot(coefficients,aes(x=year,y=beta0.5_intercept,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=beta0.025_intercept, ymax=beta0.975_intercept))+facet_wrap(~year_specific)
dev.off()

#-----------------------------------------------------------------------------------------------------------#

randomeff=unique(res[spatial==TRUE ,.(year,year_specific,quant0.5_spread,quant0.5_beta1_spread,quant0.5_beta2_spread,quant0.5_tail,quant0.5_beta3_tail,quant0.5_beta4_tail,quant0.5_range,quant0.5_sig,
                                      quant0.975_spread,quant0.975_beta1_spread,quant0.975_beta2_spread,quant0.975_tail,quant0.975_beta3_tail,quant0.975_beta4_tail,quant0.975_range,quant0.975_sig,
                                      quant0.025_spread,quant0.025_beta1_spread,quant0.025_beta2_spread,quant0.025_tail,quant0.025_beta3_tail,quant0.025_beta4_tail,quant0.025_range,quant0.025_sig)])

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//spatrange.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_range,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_range, ymax=quant0.975_range))+facet_wrap(~year_specific)+ylim(c(0,100))
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//spatsig.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_sig,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_sig, ymax=quant0.975_sig))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//spread_intercept.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_spread,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_spread, ymax=quant0.975_spread))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//spread_lon.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_beta1_spread,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_beta1_spread, ymax=quant0.975_beta1_spread))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//spread_lat.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_beta2_spread,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_beta2_spread, ymax=quant0.975_beta2_spread))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//tail_intercept.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_tail,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_tail, ymax=quant0.975_tail))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//tail_lon.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_beta3_tail,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_beta3_tail, ymax=quant0.975_beta3_tail))+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//tail_lat.pdf",width = 8, height = 4)
ggplot(randomeff,aes(x=year,y=quant0.5_beta4_tail,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+
  geom_errorbar(aes(ymin=quant0.025_beta4_tail, ymax=quant0.975_beta4_tail))+facet_wrap(~year_specific)
dev.off()

#-------------------------------------------------------------------------------------------------------------#

gevparams=unique(res[spatial==TRUE ,.(year,year_specific,mu,sigma,xi,lon,lat)])
oldgevparams=unique(res[spatial==TRUE ,.(year,year_specific,loc_gev,scale_gev,shape_gev,xi,lon,lat)])



pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//mu.pdf",width = 8, height = 4)
ggplot(gevparams,aes(x=year,y=mu,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//sig.pdf",width = 8, height = 4)
ggplot(gevparams,aes(x=year,y=sigma,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+facet_wrap(~year_specific)
dev.off()

pdf(file = "/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/params//xi.pdf",width = 8, height = 4)
ggplot(gevparams,aes(x=year,y=xi,col=year_specific))+geom_point()+geom_hline(yintercept=0,lty=2)+facet_wrap(~year_specific)
dev.off()

#-----------------------------------------------------------------------------------------------------------#

uniquelocs=unique(allres[,.(lon,lat)])

for(k in 1:(dim(uniquelocs)[1])){
  bgev_year=allres[lon== uniquelocs$lon[k] & lat==uniquelocs$lat[k] & type=="bgev_year"]
  bgev_spat=allres[lon== uniquelocs$lon[k] & lat==uniquelocs$lat[k] & type=="bgev_base"]
  normgev=allres[lon== uniquelocs$lon[k] & lat==uniquelocs$lat[k] & type=="gev"]
  
  years=unique(bgev_year$year)
  
  thisdata=amax_data[lon==uniquelocs$lon[k] & lat==uniquelocs$lat[k],]
  thisdata[,gposition:=lapply(.SD,function(x) 1/(1-gringorten(x))), .SDcols = "y"]
  
  pdf(file = paste0("/nr/project/stat/Impetus4Change/Dev_Thea/Figs/SimpleTailSpread/locyear/locnum",k,".pdf"),width = 12, height = 10)
  par(mfrow=c(3,3))
  for(yy in seq_along(years)){
    plot(bgev_year[year==years[yy],.(return_periods,return_level)],type="o",pch=21,bg="skyblue",xlim=c(0,100),ylim=c(min(allres[lon== uniquelocs$lon[k] & lat==uniquelocs$lat[k],return_level])-5,max(allres[lon== uniquelocs$lon[k] & lat==uniquelocs$lat[k] & type!= "gev",return_level])));grid()
    points(thisdata[,.(gposition,y)],pch=3)
    lines(bgev_spat[year==years[yy],.(return_periods,return_level)],type="o",pch=21,bg="yellow")
    lines(bgev_year[year==years[yy],.(return_periods,return_level)],type="o",pch=21,bg="skyblue");grid()
    
    #lines(normgev[year==years[yy],.(return_periods,return_level)],type="o",pch=21,bg="orange")
    title(paste0(years[yy]," (",round(uniquelocs$lon[k],1),", ",round(uniquelocs$lat[k],1),")"))
    points(thisdata[year==years[yy],.(gposition,y)],pch=21,bg="red")
    if(yy==1){
      legend("bottomright",col=c("skyblue","yellow","red"),c("year specific","spatial","obs"),pch=19)
    }
    
  }
  dev.off()
  
  
}

