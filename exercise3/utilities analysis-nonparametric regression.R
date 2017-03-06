source('~/GitHub/SDS383D-course-work/exercise3/gaussian process.R')
utilities <- read.csv("~/GitHub/SDS383D-course-work/exercise3/utilities.csv")

bill=as.matrix(utilities)[,2]
period=as.matrix(utilities)[,3]
temp=as.matrix(utilities)[,1]
ydata=t(as.matrix(bill/period))
xdata=t(temp)
plot(ydata~xdata)

xnew=t(as.matrix(seq(from=10,to=80)))
yhat=xnew
C=covmatrix(1:80,materncov,1,0.5,0.001)
C=solve(C)
for(i in 1:ncol(xnew)){
  Cstar=C[xnew[i],xnew[i]]
  yhat[i]=0
  for(j in 1:ncol(xdata)){
    yhat[i]=yhat[i]-C[xnew[i],xdata[j]]*ydata[j]
  }
  yhat[i]=yhat[i]/Cstar
}
lines(yhat~xnew)
