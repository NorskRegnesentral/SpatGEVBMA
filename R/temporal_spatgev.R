#' @export
make_temporal_spatgev_data=function(amax_data,to_standardize=TRUE,add_intercept=TRUE){
  Y=as.list(amax_data[,y])
  X=copy(amax_data)
  X[,stid:=NULL]
  X[,year:=NULL]
  X[,y:=NULL]
  
  S=data.table(dim1=amax_data[,year])
  
  if(to_standardize==TRUE){
    X=scale(X)
    S=scale(S)
  }
  
  X=as.matrix(X)
  S=cbind(S,0)
  
  if(add_intercept==TRUE){
    X=cbind(1,X)
    }
  
  return(list(Y=Y,X=X,S=S))
}
