source('~/GitHub/SDS383D-course-work/exercise3/local polynomial regression/local linear estimator.R')
utilities <- read.csv("~/GitHub/SDS383D-course-work/exercise3/utilities.csv")

bill=as.matrix(utilities)[,2]
period=as.matrix(utilities)[,3]
temp=as.matrix(utilities)[,1]
ydata=t(as.matrix(bill/period))
xdata=t(temp)
plot(ydata~xdata,xlab="temprature",ylab="average bill")

xgrid=t(as.matrix(seq(10,80,by=0.1)))
yhat=fitsmoother(xgrid,ydata,xdata,3)
lines(xgrid,yhat,col="red")

width=seq(3,15,by=0.1)
loocverror=width
for(i in 1:length(width)){
  loocverror[i]=loocv(ydata,xdata,width[i])
}
plot(loocverror~width,type="l")

yhat=fitsmoother(xgrid,ydata,xdata,6.8)
plot(ydata~xdata,xlab="temprature",ylab="average bill")
lines(xgrid,yhat,col="red")

ydatahat=fitsmoother(xdata,ydata,xdata,6.8)
logresidue=2*log(abs(ydata-ydatahat))

loocverror=width
for(i in 1:length(width)){
  loocverror[i]=loocv(logresidue,xdata,width[i])
}
plot(loocverror~width,type="l")

plot(logresidue~xdata,xlab="temprature",ylab="logerror")
lrhat=fitsmoother(xgrid,logresidue,xdata,9)
lines(xgrid,lrhat,col="red")

dy=sqrt(exp(lrhat))*1.96

plot(ydata~xdata,xlab="temprature",ylab="average bill")
lines(xgrid,yhat,col="red")
lines(yhat+dy~xgrid,col="blue")
lines(yhat-dy~xgrid,col="blue")