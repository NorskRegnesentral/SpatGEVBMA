#' @param amax_data 
#'
#' @param to_standardize 
#' @param add_intercept 
#'
#' @export
make_temporal_spatgev_data=function(amax_data,to_standardize=TRUE,add_intercept=TRUE,tau_effect="temporal"){
  #tau_effect er "temporal" eller "spatial".
  
  Y=as.list(amax_data[,y])
  X=copy(amax_data)
  X[,stid:=NULL]
  X[,year:=NULL]
  X[,y:=NULL]
  
  if(tau_effect=="temporal"){
    S=data.table(dim1=amax_data[,year])
  }
  
  if(tau_effect=="spatial"){
    S=data.table(dim1=amax_data[,lon],dim2=amax_data[,lat])
    
  }
  
  
  
  if(to_standardize==TRUE){
    X=scale(X)
    S=scale(S)
  }
  
  X=as.matrix(X)
  
  if(tau_effect!="spatial"){
    S=cbind(S,0)
  }
  
  if(add_intercept==TRUE){
    X=cbind(1,X)
    }
  
  return(list(Y=Y,X=X,S=S))
}
