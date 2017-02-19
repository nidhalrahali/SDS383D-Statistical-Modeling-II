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

