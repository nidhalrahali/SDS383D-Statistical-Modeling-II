source('~/GitHub/SDS383D-course-work/exercise3/local linear estimator.R')
utilities <- read.csv("~/GitHub/SDS383D-course-work/exercise3/utilities.csv")

bill=as.matrix(utilities)[,2]
period=as.matrix(utilities)[,3]
temp=as.matrix(utilities)[,1]
y=t(as.matrix(bill/period))
x=t(temp)
plot(y~x)

xnew=t(as.matrix(seq(10,80,by=0.1)))
yhat=fitsmoother(xnew,y,x,3)
lines(xnew,yhat,col="red")

width=seq(3,15,by=0.1)
loocverror=width
for(i in 1:length(width)){
  loocverror[i]=loocv(y,x,width[i])
}
plot(loocverror~width,type="l")

xnew=t(as.matrix(seq(10,80,by=0.1)))
yhat=fitsmoother(xnew,y,x,7)
lines(xnew,yhat,col="red")

yhat=fitsmoother(x,y,x,7)
residue=y-yhat