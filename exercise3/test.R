n=1000

x=matrix(ncol=n)
for(i in 1:n){
  x[i]=-5+10/n*i
}
xtrain=sample(x,400,replace=TRUE)
xtest=sample(x,600,replace=TRUE)

f=function(x){
  return=x^3/10-x
}
ytrain=f(xtrain)
ytest=f(xtest)

h=matrix(ncol=100)
for(i in 1:100){
  h[i]=i/50
}
sqerror=matrix(ncol=100)
for(i in 1:100){
error=testmodel(ytrain,xtrain,ytest,xtest,h[i],"gaussian")
sqerror[i]=tcrossprod(error,y=NULL)/ncol(error)
}
plot(sqerror~h)
