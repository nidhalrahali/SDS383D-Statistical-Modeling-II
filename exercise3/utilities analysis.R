source('~/GitHub/SDS383D-course-work/exercise3/local linear estimator.R')
utilities <- read.csv("~/GitHub/SDS383D-course-work/exercise3/utilities.csv")

bill=as.matrix(utilities)[,2]
period=as.matrix(utilities)[,3]
temp=as.matrix(utilities)[,1]
ydata=t(as.matrix(bill/period))
xdata=t(temp)
plot(ydata~xdata)

xnew=t(as.matrix(seq(10,80,by=0.1)))
yhat=fitsmoother(xnew,ydata,xdata,3)
lines(xnew,yhat,col="red")

width=seq(3,15,by=0.1)
loocverror=width
for(i in 1:length(width)){
  loocverror[i]=loocv(ydata,xdata,width[i])
}
plot(loocverror~width,type="l")

xnew=t(as.matrix(seq(10,80,by=0.1)))
yhat=fitsmoother(xnew,ydata,xdata,6.8)
plot(ydata~xdata)
lines(xnew,yhat,col="red")

ydatahat=fitsmoother(xdata,ydata,xdata,6.8)
#ydataloghat=fitsmoother(xdata,log(ydata),xdata,6.8)
residue=ydata-ydatahat
plot(residue~xdata)
#logresidue=log(ydata)-ydataloghat
#plot(logresidue~xdata)
H=smoothingmatrix(xdata,6.8)
sigma=sd(logresidue)*(length(xdata)-1)/(length(xdata)-2*sum(diag(H))+sum(diag(crossprod(H,H))))

deltay=xnew
for(i in 1:ncol(xnew)){
  w=weight(xnew[i],xdata,6.8)
  w=w/sum(w)
  deltay[i]=1.96*sqrt(tcrossprod(w,w))*sigma
}
plot(ydata~xdata)
lines(yhat+deltay~xnew)
lines(yhat-deltay~xnew)