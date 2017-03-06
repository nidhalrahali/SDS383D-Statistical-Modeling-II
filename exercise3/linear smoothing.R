# this script generate a nonlinear function with gaussian noise
# then use linear smoother to recover the function

source('~/GitHub/SDS383D-course-work/exercise3/weight.R')

# the function we use to generate data
f=function(x){
  return=x^3/10-x
}

n=500
x=matrix(nrow=1,ncol=n)
for(i in 1:n){
  x[1,i]=-5+10/n*i
}

# curve of original function
y=f(x)
plot(y~x,type="l")

# data points after noise
ynoise=f(x)+rnorm(n,sd=0.5)
points(x,ynoise)

# normalize y, x was already normalized manually

# use linear smoother to fit the data, (y,x) input data, h: bandwidth
# ty:type of kernel, cl: color used to draw a line

fitsmoother=function(y,x,h,ty,cl){
  yhat=matrix(ncol=n)
  for(i in 1:ncol(y)){
    yhat[i]=y %*% t(weight(x[i],x,h,type=ty))
  }
  lines(x,yhat,col=cl)
}

# compare different bandwidth and kernel
fitsmoother(ynoise,x,100,"gaussian","red")
fitsmoother(ynoise,x,10,"gaussian","green")
fitsmoother(ynoise,x,10,"uniform","blue")
fitsmoother(ynoise,x,5,"gaussian","purple")

