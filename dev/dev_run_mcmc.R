rm(list = ls())

library(SpatGEVBMA)

output.path <- "~/Desktop/toAlex/"
covariates.folder <- "~/Desktop/toAlex/Data/Cov/"
station.annualMax.file <- "~/Desktop/toAlex/Data/Obs/AM_final.xlsx"
return.period <- c(2,5,10,20,25,50,100,200)
post.quantiles <- c(0.025,0.5,0.975)
show.uncertainty <- TRUE
coordinate.type <- "XY"
table.format = "html"
mcmc.reps <- 2e2 
burn.in <- 2e1
cores <- 16 # 20
annualMax.name <- NULL
create.tempfiles <- TRUE
keep.temp.files <- FALSE
save.all.output <- TRUE
testing <- FALSE
transform.output = "UTM_33_to_LatLon"
fixed.xi = 0

station.annualMax.sheet <- 12
station.locations.file <- "~/Desktop/toAlex/Data/Meta/metadata_stations_plu_and_geo.txt"
output.folder.name <- "SpatGEV.res.1440min"



print_map = function(Q,
                     shortName,
                     longName,
                     post.quantiles,
                     IQRLongName,
                     show.uncertainty = TRUE,
                     ww.na = NULL,
                     n = dim(Q)[1],
                     output.name,
                     nx,
                     ny)
{
  ## Default variables attributes
  output.unit = "millimeter"  
  output.missval <- -999.99
  output.chunksizes <- c(length(output.x)/3,length(output.y)/3)

  ## Default global attributes  
  output.references <- "Output from the the SpatGEVBMA.wrapper function in the R-package SpatGEVBMA, developed in Dyrrdal, A. V. et al. (2015)"
  output.referencesRpackage <- "https://github.com/NorskRegnesentral/SpatGEVBMA"
  output.referencesPaper <- "Dyrrdal, A. V., Lenkoski, A., Thorarinsdottir, T. L., & Stordal, F. (2015). Bayesian hierarchical modeling of extreme hourly precipitation in Norway. Environmetrics, 26(2), 89-106."
  output.var.version <- "no.version"  
                                        #output.proj4.string <- "+proj=utm+zone=33+ellps=WGS84"  # MJ: This shouldn't really be hardcoded...
  output.prod.date = substr(Sys.time(), 1, 10)
                                        #output.conventions = "CF-1.4"
  output.institution = "Norwegian Computing Center (Norsk Regnesentral)"

  ncvar_defList <- list()

  for (i in 1:length(longName))
  {
    ncvar_defList[[i]] <- ncvar_def(name=shortName[i],longname=longName[i],
                                    units=output.unit,
                                    dim=dim.list,
                                    missval=output.missval, # How missing values in input data are defined 
                                    chunksizes = output.chunksizes)
  }
      

  if (show.uncertainty)
  {
    this.IQR <- length(post.quantiles)+1
    shortName <- c(shortName,"IQR")
    longName <- c(longName,IQRLongName)
    ncvar_defList[[this.IQR]] <- ncvar_def(name=shortName[this.IQR],longname=longName[this.IQR],
                                           units=output.unit,
                                           dim=dim.list,
                                           missval=output.missval, # How missing values in input data are defined 
                                           chunksizes = output.chunksizes)
  }

  notNA <- which(!(1:n %in% ww.na))

  full.Z.p <- matrix(NA,ncol=n,nrow=dim(Z)[2])
  full.Z.p[,notNA] <- Q

        retMat.quant=list()
      for (r in 1:length(post.quantiles))
        {
        retMat.quant[[r]] <- matrix(full.Z.p[r,], ncol=ny,nrow=nx)
        if (!is.null(transform.output)){
          original.image$z <- retMat.quant[[r]]
          interp.values <- interp.surface(obj=original.image, loc=XYGrid)
          retMat.quant[[r]] <- matrix(interp.values, ncol=ny,nrow=nx,byrow=T)
        }
          ncvar_put(outputNc,varid=ncvar_defList[[r]],vals=c(retMat.quant[[r]][,ny:1]))  ##check
        }
      if (show.uncertainty)
        {
          retMat.IQR <- matrix(full.Z.p[this.IQR+1,]-full.Z.p[this.IQR,],ncol=ny,nrow=nx) # This part of full.Z.p is not transformed above
          if (!is.null(transform.output)){
            original.image$z <- retMat.IQR
            interp.values <- interp.surface(obj=original.image, loc=XYGrid)
            retMat.IQR <- matrix(interp.values, ncol=ny,nrow=nx,byrow=T)
          }
          ncvar_put(outputNc,varid=ncvar_defList[[this.IQR]],vals=c(retMat.IQR[,ny:1]))
        }

  ## Putting global attributes to the nc-file:
  ## Default global attributes  
  ncatt_put(outputNc,varid=0, attname = "references", attval = output.references)
  ncatt_put(outputNc,varid=0, attname = "Rpackage", attval = output.referencesRpackage)
  ncatt_put(outputNc,varid=0, attname = "Paper", attval = output.referencesPaper)
  ncatt_put(outputNc,varid=0, attname = "var.version", attval = output.var.version)
  ##ncatt_put(outputNc,varid=0, attname = "proj4.string", attval = output.proj4.string)
  ncatt_put(outputNc,varid=0, attname = "prod.date", attval = output.prod.date)
  ##ncatt_put(outputNc,varid=0, attname = "conventions", attval = output.conventions)
  ncatt_put(outputNc,varid=0, attname = "institution", attval = output.institution)
  
  nc_close(outputNc)

        if (coordinate.type=="XY"){
        if (is.null(transform.output)){
          lab.name=c("x-coord","y-coord")
        } else {
          lab.name=c("Lon","Lat")
        }
      }
      if (coordinate.type=="LatLon"){
        lab.name=c("Lon","Lat")
      }
      filename.pdf <- file.path(output.folder,paste("posterior.return.level.",return.period[j],"grid.pdf",sep=""))
      pdf(filename.pdf,width=7, height=7)
      for (r in 1:length(post.quantiles))
        {
          image.plot(output.x,output.y,retMat.quant[[r]],main=paste("Posterior ", post.quantiles[r], "-quantile \n ", return.period[j]," year return value with ", annualMax.name," data",sep=""),xlab=lab.name[1],ylab=lab.name[2])
          points(S[,1],S[,2])
        }
      if (show.uncertainty)
        {
          image.plot(output.x,output.y,retMat.IQR,main=paste("Interquartile range uncertainty plot \n ", return.period[j]," year return value with ", annualMax.name," data",sep=""),xlab=lab.name[1],ylab=lab.name[2])
          points(S[,1],S[,2])
        }
      dev.off()

}
