library(mvtnorm)
library(LaplacesDemon)
gibbssampler=function(xx,xy,yy,p,t){
  n=length(yy)
  invsig2=rep(1,n)
  m=rep(0,4)
  mu=m
  beta=matrix(0,ncol=4,nrow=n)
  lambdastar=1+n
  nustar=2*n
  sigma=diag(1,nrow=4,ncol=4)
  psi=diag(1,nrow=4,ncol=4)
  for(i in 1:500){
    sigmainv=solve(sigma)
    temp=matrix(0,nrow=4,ncol=4)
    for(i in 1:n){
      temp=temp+xx[i,,]*invsig2[i]
    }
    sigmastar=solve(temp+sigmainv)
    temp=matrix(0,nrow=4,ncol=1)
    for(i in 1:n){
      temp=temp+t(xy)[i,]*invsig2[i]
    }
    mustar=sigmastar%*%(temp+sigmainv%*%mu)
    betamean=colMeans(beta)
    mstar=(mu+betamean*n)/lambdastar
    temp=matrix(0,nrow=4,ncol=4)
    psistar=psi+crossprod(beta-betamean,beta-betamean)+n*crossprod(beta-mu,beta-mu)/lambdastar
    for(i in 1:n){
      r=(yy[i]-2*crossprod(beta[i,],xy[,i])+crossprod(beta[i,],xx[i,,])%*%beta[i,])/2
      invsig2[i]=rgamma(1,shape=p[i]/2+1,rate=r)
    }
    beta=rmvnorm(n,mustar,sigmastar)
    mu=as.vector(rmvnorm(1,m,sigma))
    sigma=rinvwishart(nustar,psistar)
  }
}