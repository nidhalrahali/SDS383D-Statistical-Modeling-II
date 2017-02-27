kernel=function(d){
  return=exp(-d*d/2)/sqrt(2*pi)
}

weight=function(x,xdata,h){
  d=xdata-x
  s1=kernel(d/h)*d
  s2=s1*d
  s1=sum(s1)
  s2=sum(s2)
  return=kernel(d/h)*(s2-s1*d)
}

fitsmoother=function(xnew,y,x,h){
  ymean=mean(y)
  ynormalize=y-ymean
  yhat=xnew
  for(i in 1:ncol(xnew)){
    w=weight(xnew[i],x,h)
    w=w/sum(w)
    yhat[i]=tcrossprod(ynormalize,w)+ymean
  }
  return=yhat
}

loocv=function(y,x,h){
  yhat=fitsmoother(x,y,x,h)
  sqerror=0
  for(i in 1:ncol(x)){
    w=weight(x[i],x,h)
    w=w[i]/sum(w)
    sqerror=sqerror+((y[i]-yhat[i])/(1-w))^2
  }
  return=sqerror
}