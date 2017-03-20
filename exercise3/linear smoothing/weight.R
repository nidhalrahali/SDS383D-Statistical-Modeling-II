# compute weight function for linear smooth
# d: the distance between data points
# h: the bandwidth
# type: type of kernel, "uniform" uniform kernel, "gaussian" for gaussian kernel,default is "uniform"

weight=function(xhat,x,h,type="uniform"){
  d=x-matrix(xhat,ncol=ncol(x))
  if(type=="gaussian"){
    r=exp(-d*d/2)/sqrt(2*pi)/h
  }
  if(type=="uniform"){
    r=matrix(nrow=1,ncol=ncol(d))
    for(i in 1:ncol(d)){
      if(abs(d[i])<=1){
        r[i]=0.5/h
      }
      else{
        r[i]=0
      }
    }
  }
  return=r
}

predictionerror=function(ytrain,xtrain,ytest,xtest,h,t){
  #normalize all data
  ymean=mean(ytrain)
  ynorm=ytrain-matrix(ymean,ncol=length(ytrain))
  xnorm=xtrain-matrix(mean(xtrain),ncol=length(xtrain))
  xtestnorm=xtest-matrix(mean(xtrain),ncol=length(xtest))
  
  # fit smoother and record the error
  error=matrix(ncol=length(ytest))
  for(i in 1:length(ytest)){
    error[i]=ynorm %*% t(weight(xtestnorm[i],xnorm,h,type=t))+ymean-ytest[i]
  }
  return(error)
}
