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

covmatrix=function(x,covfun,b,t1,t2){
  n=length(x)
  covmatrix=matrix(ncol=n,nrow=n)
  for(i in 1:n){
    for (j in i:n){
      covmatrix[i,j]=covfun(abs(x[i]-x[j]),b,t1,t2)
      covmatrix[j,i]=covmatrix[i,j]
    }
  }
  return=covmatrix
}

fitsmoother=function(xnew,y,x,covfun,b,t1,t2){
  yhat=xnew
  for(i in 1:ncol(xnew)){
    C=solve(covmatrix(xnew[i],x,covfun,b,t1,t2))
    w=C[1,2:ncol(C)]/C[1,1]
    yhat[i]=-tcrossprod(y,w)
  }
  return=yhat
}