rcp45=nc_open("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/inference/rcp45/1991/MAP.nc")
tas1=ncvar_get(rcp45,"MAP")

image.plot(tas1)


rcp45=nc_open("/nr/project/stat/ClimDesign/WP3/Data/fromOskar/prediction/rcp45/2071/MAP.nc")
tas2=ncvar_get(rcp45,"MAP")

diff=tas2-tas1

diff=diff[is.na(diff)==FALSE]

hist(diff[1:10000])
