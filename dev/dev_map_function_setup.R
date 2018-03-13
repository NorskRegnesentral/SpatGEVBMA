rm(list = ls())

library(SpatGEVBMA)

setwd("~/Desktop/toAlex/SpatGEV.res.1440min")

load("./Temp/temp_checkpoint_0.RData")
load("./Temp/temp_checkpoint_1.RData")
load("./Temp/temp_checkpoint_2.RData")
load("./Temp/temp_checkpoint_3.RData")
load("./Temp/temp_checkpoint_4.RData")

                                        # Checkpoint 4

## Transforming coordinate system output if applicable

## Need one of the input files to specify parameters in the ncdf output file
ncnc <- nc_open(cov.files.path[1])


if(!is.null(transform.output))
{
  UTM.zone <- as.numeric(substr(transform.output,start=5,stop=6))
  
                                        # Assume coordinate.type=="XY"
  rangeX <- range(indX)
  rangeY <- range(indY)
  
  allXYMat <- data.frame(X=allX,Y=allY)
  
                                        #Set projection and zone
  attr(allXYMat, "projection") <- "UTM"
  attr(allXYMat, "zone") <- UTM.zone
  
                                        #Compute LL coordinates
  allLLMat <- as.matrix(round(convUL(allXYMat, km=FALSE), digits=4))
  
  indLon <- seq(from = min(allLLMat[,1]),to = max(allLLMat[,1]),length.out = nx)
  indLat <- seq(from = min(allLLMat[,2]),to = max(allLLMat[,2]),length.out = ny)
  
  allLon <- rep(indLon,each=ny)
  allLat <- rep(indLat,times=nx)
  
  allLL <- cbind(X=allLon,Y=allLat)
    
    # Set projection and zone
  attr(allLL, "projection") <- "LL"
  attr(allLL, "zone") <- UTM.zone
  
                                        #Compute UTM coordinates
  XYGrid <- as.matrix(round(convUL(allLL, km=FALSE), digits=0))
  
                                        # Transform the S matrix to LatLon as well
  
  S.new <- cbind(X=S[,1],Y=S[,2])
  
  attr(S.new, "projection") <- "UTM"
  attr(S.new, "zone") <- UTM.zone
  
                                        #Compute LL coordinates
  S <- as.matrix(round(convUL(S.new, km=FALSE), digits=4))
  
}

## Writing final results to netcdf-file and image plots

## Define the dimensions
if (coordinate.type=="XY")
{
  if (is.null(transform.output))
  {
    
    output.x <- indX
    output.y <- indY
    
    x.ncdf <- ncdim_def( "X", "meters", output.x)
    y.ncdf <- ncdim_def( "Y", "meters", output.y[ny:1])
    
    dim.list <- list(X=x.ncdf,Y=y.ncdf)
  } else {
    
    original.image <- list(x=indX,y=indY)
    
    output.x <- indLon
    output.y <- indLat
    
    x.ncdf <- ncdim_def( "Lon", "degrees_E", output.x)
    y.ncdf <- ncdim_def( "Lat", "degrees_N", output.y[ny:1])
    
    dim.list <- list(Lon=x.ncdf,Lat=y.ncdf)
    
  }
  
}

if (coordinate.type=="LatLon")
{
  output.x <- indX
  output.y <- indY
  
  x.ncdf <- ncdim_def( "Lon", "degrees_E", output.x)
  y.ncdf <- ncdim_def( "Lat", "degrees_N", output.y[ny:1])
  dim.list <- list(Lon=x.ncdf,Lat=y.ncdf)
}


j = 1
## Just general definitions
shortName <- paste("quant_",gsub(".","_",post.quantiles,fixed=TRUE),sep="")  # The gsub thing replaces the dot with a underscore
longName <- paste(post.quantiles," quantile of the marginal ",
                  "posterior distribution for the maximum precipition over ",
                  return.period[j]," years based on data: ",
                  annualMax.name,".",sep="")
w_median = which(post.quantiles == 0.5)
if(length(w_median) > 0)
{
  longName[w_median] <- paste("Median of the marginal posterior ",
                              "distribution for the ",return.period[j],
                              " return level for precipitation based on data: ",
                              annualMax.name,".",sep="")
}

IQRLongName <- paste("Interquartile range uncertainty measure: Difference" ,
                     "between 0.75-quantile and 0.25-quantile for ",
                     "the maximum precipitaion over ",
                     return.period[j], " years based on data: ",
                     annualMax.name,".",sep="")

filename.nc <- paste0(output.folder,"/posterior.grid_return_",return.period[j],".nc",sep="")
filename.pdf <- paste0(output.folder,paste("/posterior.return.level.",return.period[j],"grid.pdf",sep=""))

main.quantile = paste("Posterior ", post.quantiles, "-quantile \n ", return.period[j]," year return value with ", annualMax.name," data",sep="")
main.iqr = paste("Interquartile range uncertainty plot \n ", return.period[j]," year return value with ", annualMax.name," data",sep="")
output.name <- paste(filename.nc,"_return_",return.period[j],".nc",sep="")

cov <- as.matrix(gridData$covariates)

ww.na <- which(apply(is.na(cov),1,"any"))

load("./imputation.RData")

save(Z.p, shortName, longName, post.quantiles,
     IQRLongName, ww.na, n, output.x, output.y, output.name,
     filename.nc, filename.pdf, main.quantile, main.iqr,
     nx,ny,dim.list,all.post.quantiles,
     transform.output,original.image, XYGrid,coordinate.type,S, file = "./Temp/print_maps_injection.RData")
