library(mvtnorm)
library(MCMCpack)
gibbssampler=function(xx,xy,yy,t,s){
  n=length(yy)
  invsig2=1
  m=rep(0,4)
  mu=m
  beta=matrix(0,ncol=4,nrow=n)
  lambda=1
  lambdastar=lambda+n
  nu=4
  nustar=n+nu
  sigma=diag(4)
  psi=diag(1,nrow=4,ncol=4)
  retbeta=array(dim=c(n,t,4))
  retsigma=array(dim=c(t,4,4))
  retsiginv=rep(0,t)
  for(j in 1:(500+t)){
    sigmainv=solve(sigma)
    for(i in 1:n){
      sigmastar=solve(sigmainv+xx[i,,]*invsig2)
      mustar=sigmastar%*%(sigmainv%*%mu+xy[,i]*invsig2)
      beta[i,]=rmvnorm(1,mustar,sigmastar)
      if(j>500){
        retbeta[i,j-500,]=beta[i,]
      }
    }
    betamean=colMeans(beta)
    mstar=(m*lambda+betamean*n)/lambdastar
    psistar=psi
    for(i in 1:n){
      psistar=psi+tcrossprod(beta[i,]-betamean,beta[i,]-betamean)+tcrossprod(beta[i,]-mu,beta[i,]-mu)*n*lambda/lambdastar
    }
    r=0
    for(i in 1:n){
      r=r+(yy[i]-2*crossprod(beta[i,],xy[,i])+crossprod(beta[i,],xx[i,,])%*%beta[i,])/2
    }      
    invsig2=rgamma(1,shape=s/2+1,rate=r)
    if(j>500)retsiginv[j-500]=invsig2
    mu=as.vector(rmvnorm(1,mstar,sigma))
    sigma=riwish(nustar,psistar)
    if(j>500)retsigma[j-500,,]=sigma
  }
  return=list(betasample=retbeta,sigmasample=retsigma,siginvsample=retsiginv)
}