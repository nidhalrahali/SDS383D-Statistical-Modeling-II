source('~/GitHub/SDS383D-course-work/exercise3/gaussian process.R')
utilities <- read.csv("~/GitHub/SDS383D-course-work/exercise3/utilities.csv")

#read data
bill=as.matrix(utilities)[,2]
period=as.matrix(utilities)[,3]
temp=as.matrix(utilities)[,1]
ydata=t(as.matrix(bill/period))
xdata=t(temp)
plot(ydata~xdata,xlab='temperature',ylab='bill')

#parameters
b=30
t1=10
t2=0.001
sigsquare=0.01


xnew=t(as.matrix(seq(from=10,to=80)))
yhat=xnew
yvariance=xnew

C=covmatrix(xdata,materncov,b,t1,t2)
D=solve(C+sigsquare*diag(ncol(C)))
cstarstar=materncov(0,b,t1,t2)

for(i in 1:ncol(xnew)){
  cstar=xdata
  for(j in 1:ncol(xdata)){
    cstar[j]=materncov(abs(xnew[i]-xdata[j]),b,t1,t2)
  }
    W=tcrossprod(cstar,D)
    yhat[i]=tcrossprod(ydata,W)
    yvariance[i]=sqrt((cstarstar-tcrossprod(cstar,W)))*1.96
}

#plot
plot(ydata~xdata,xlab='temperature',ylab='bill')
lines(yhat~xnew,col='red')
lines(yhat+yvariance~xnew,col='blue')
lines(yhat-yvariance~xnew,col='blue')

sigsquare=0.01
t2=0
t1range=seq(1,100,by=1)
brange=seq(1,100,by=1)
likelihood=matrix(ncol=length(brange),nrow=length(t1range))
for(i in 1:length(t1range)){
  for(j in 1:length(brange)){
    C=covmatrix(xdata,materncov,brange[j],t1range[i],t2)
    C=C+diag(nrow(C))*sigsquare
    likelihood[i,j]=logml(ydata,C)
  }
}
contour(t1range,brange,likelihood)
