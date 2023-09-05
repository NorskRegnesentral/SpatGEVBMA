library(ncdf4)
library(fields)
library(omnibus)

a=nc_open("/nr/project/stat/ClimDesign/WP3/Res/rcp45/prediction/2071/res_180min_rcp45/posterior.grid_return_50.nc")
quantis=mirror(ncvar_get(a,"quant_0_5"))
nc_close(a)


a=nc_open("/nr/project/stat/ClimDesign/WP3/Res/rcp45/inference/res_180min_rcp45/posterior.grid_return_50.nc")
quantis2=mirror(ncvar_get(a,"quant_0_5"))
nc_close(a)

image.plot(quantis)#,zlim=c(10,55))
image.plot(quantis2)#,zlim=c(10,55))

image.plot(quantis2-quantis,zlim=c(-2,2))

#-------------------------------------------------------------------------

a=nc_open("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/inference/rcp45/1991/MSPJJA.nc")
var1=ncvar_get(a,"MSPJJA")

a=nc_open("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/prediction/rcp45/2041/MSPJJA.nc")
var2=ncvar_get(a,"MSPJJA")

image.plot(t(flip(t(flip(flip(var2-var1))))))



a=nc_open("/nr/project/stat/ClimDesign/WP3/Res/rcp45/prediction/2071/res_10min_rcp45/posterior.grid_return_2.nc")
quantis=ncvar_get(a,"quant_0_5")

image.plot(quantis)


#------------------------------------------------#

a=nc_open("/nr/project/stat/ClimDesign/WP3/Res/rcp45/inference/res_10min_rcp45/posterior.grid_return_2.nc")
quantis2=t(flip(t(flip(flip(ncvar_get(a,"quant_0_5"))))))
nc_close(a)

image.plot(quantis2)

image.plot(t(flip(t(flip(flip(quantis-quantis2))))))

