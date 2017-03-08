source('~/GitHub/SDS383D-course-work/exercise3/gaussian process.R')
utilities <- read.csv("~/GitHub/SDS383D-course-work/exercise3/utilities.csv")

bill=as.matrix(utilities)[,2]
period=as.matrix(utilities)[,3]
temp=as.matrix(utilities)[,1]
ydata=t(as.matrix(bill/period))
xdata=t(temp)
plot(ydata~xdata)

b=3
t1=0.5
t2=0.001
xnew=t(as.matrix(seq(from=10,to=80)))
yhat=xnew
C=covmatrix(xdata,materncov,b,t1,t2)
sigsquare=0.01
D=solve(C+diag(x=sigsquare,nrow=nrow(C),ncol=ncol(C)))
for(i in 1:ncol(xnew)){
  cstar=xdata
  for(j in 1:ncol(xdata)){
    cstar[j]=materncov(abs(xnew[i]-xdata[j]),b,t1,t2)
  }
    W=tcrossprod(cstar,D)
    yhat[i]=tcrossprod(ydata,W)
}
plot(ydata~xdata)
lines(yhat~xnew,col='red')


t2=0
t1range=seq(100,200,by=10)
brange=seq(0.5,4,by=0.5)
likelihood=matrix(ncol=length(brange),nrow=length(t1range))
for(i in 1:length(t1range)){
  for(j in 1:length(brange)){
    C=covmatrix(xdata,materncov,brange[j],t1range[i],t2)
    D=solve(C+diag(x=sigsquare,nrow=nrow(C),ncol=ncol(C)))
    likelihood[i,j]=logml(ydata,D)
  }
}
C=covmatrix(xdata,materncov,60,6,t2)
D=solve(C+diag(x=sigsquare,nrow=nrow(C),ncol=ncol(C)))
a=logml(ydata,D)