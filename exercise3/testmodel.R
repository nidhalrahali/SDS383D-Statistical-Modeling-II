source('~/GitHub/SDS383D-course-work/exercise3/predictionerror.R')

# smooth
f=function(x){
  return=x^3/10-x
}
plot(f)

# wiggly
g=function(x){
  return=sin(x*100)+x
}
plot(g)

# number of samples
n=100

x=matrix(ncol=n)
for(i in 1:n){
  x[i]=-5+10/n*i
}

# seperate into training and test
xrand=sample(x,n)
xtrain=xrand[1:(n*2/3)]
xtest=xrand[(n*2/3+1):n]

# smooth and quiet
ytrain=f(xtrain)+rnorm(length(xtrain),sd=0.5)
ytest=f(xtest)+rnorm(length(xtest),sd=0.5)

# smooth and noisy
ytrain=f(xtrain)+rnorm(length(xtrain),sd=2)
ytest=f(xtest)+rnorm(length(xtest),sd=2)

# wiggly and quiet
ytrain=g(xtrain)+rnorm(length(xtrain),sd=0.5)
ytest=g(xtest)+rnorm(length(xtest),sd=0.5)

# wiggly and noisy
ytrain=g(xtrain)+rnorm(length(xtrain),sd=2)
ytest=g(xtest)+rnorm(length(xtest),sd=2)


plot(ytrain~xtrain)

h=matrix(ncol=100)
for(i in 1:100){
  h[i]=i/10+2
}
sqerror=matrix(ncol=100)
for(i in 1:100){
error=predictionerror(ytrain,xtrain,ytest,xtest,h[i],"gaussian")
sqerror[i]=tcrossprod(error,y=NULL)/ncol(error)
}
plot(sqerror~h)
