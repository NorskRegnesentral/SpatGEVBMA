##########################################################
# 
# Converts lon and lat coordinates to UTM-33N
#
# LLfile: txt file with LL coordinates and header X Y
#         format: "path to file"
# outpath: "path to where results should be written"
#
# Anita Verpe Dyrrdal, met.no, 2009
#
##########################################################

LLtoUTM <- function(LLfile, outpath) {

library(PBSmapping)

#Read file with lon/lat coordinates
if(is.list(LLfile)) {
	LL <- t(sapply(LLfile, unlist))
} else {
	LL <- as.matrix(read.table(LLfile, header = TRUE))
}

#Number of locations
no.loc <- nrow(LL)

#Set projection and zone
attr(LL, "projection") <- "LL"
attr(LL, "zone") <- 33

#Compute UTM coordinates
UTM <- as.matrix(round(convUL(LL, km=FALSE), digits=0))

#Write out UTM file
write.table(UTM, paste(outpath, "coordUTM.txt",sep=""), quote=FALSE, row.names=FALSE)

}

UTMtoLL <- function(UTMfile, outpath) {
  
  library(PBSmapping)
  
  #Read file with UTM coordinates
  UTM <- as.matrix(read.table(UTMfile, header = TRUE))
  
  #Number of locations
  no.loc <- nrow(UTM)
  
  #Set projection and zone
  attr(UTM, "projection") <- "UTM"
  attr(UTM, "zone") <- 33
  
  #Compute LL coordinates
  LL <- as.matrix(round(convUL(UTM, km=FALSE), digits=4))
  
 # return(LL)
  
  #Write out LL file
  write.table(LL, paste(outpath, "coordLL.txt",sep=""), quote=FALSE, row.names=FALSE)
  
}











