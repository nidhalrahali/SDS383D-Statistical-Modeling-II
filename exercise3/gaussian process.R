library(mvtnorm)

materncov=function(d,b,t1,t2){
  r=t1^2*(1+sqrt(5)*d/b+5*d^2/(3*b^2))*exp(-sqrt(5)*d/b)
  if(d==0){
    r=r+t2^2
  }
  return=r
}

secov=function(d,b,t1,t2){
  r=t1^2**exp(-d^2/(2*b^2))
  if(d==0){
    r=r+t2^2
  }
  return=r
}

gaussian=function(covfun,x,b,t1,t2){
  covariance=matrix(ncol=length(x),nrow=length(x))
  for(i in 1:length(x)){
    for(j in i:length(x)){
      covariance[i,j]=covfun(x[j]-x[i],b,t1,t2)
      covariance[j,i]=covariance[i,j]
    }
  }
  return=rmvnorm(1,rep(0,length(x)),covariance)
}