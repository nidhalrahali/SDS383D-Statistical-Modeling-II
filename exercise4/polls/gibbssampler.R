library(mvtnorm)
library(truncnorm)

gibbssampler=function(y,x,state,beta,mu,g,lambdasq,sigsq,tausq,d,eta,t){
  n=length(y)
  f=length(beta)
  m=length(mu)
  z=rep(0,n)
  sigsqsample=rep(0,t)
  muzero=0
  beta=t(as.matrix(beta))
  betasample=matrix(nrow=t,ncol=f)
  betavarsample=array(dim=c(t,f,f))
  musample=matrix(nrow=t,ncol=m)
  muvarsample=matrix(nrow=t,ncol=m)
  for(i in 1:(1000+t)){
    if(i%%100==0){
      print(i)
    }
    for(j in 1:n){
      zmean=tcrossprod(x[j,], beta)+mu[state[j]]
      if(y[j]==1){z[j]=rtruncnorm(1,a=0,mean=zmean,sd=sqrt(sigsq))}
      if(y[j]==0){z[j]=rtruncnorm(1,b=0,mean=zmean,sd=sqrt(sigsq))}
    } 
    zmu=z
    for(j in 1:n){
      zmu[j]=z[j]-mu[state[j]]
    }
    betavar=diag(1,f)/tausq+crossprod(x,x)/sigsq
    betavar=solve(betavar)
    betamean=betavar%*%crossprod(x,zmu)/sigsq
    beta=rmvnorm(1,betamean,betavar)
    if(i>1000){
      betasample[i-1000,]=beta
      betavarsample[i-1000,,]=betavar
    }
    mumean=rep(0,m)
    for(j in 1:n){
      mumean[state[j]]=mumean[state[j]]+z[j]-tcrossprod(x[j,], beta)
    }
    mumean=(mumean*lambdasq+muzero*sigsq)/(sigsq+lambdasq*g)
    muvar=lambdasq*sigsq/(sigsq+lambdasq*g)
    mu=rnorm(m,mumean,sqrt(muvar))
    if(i>1000){
      musample[i-1000,]=mu
      muvarsample[i-1000,]=muvar
    }
    betarms=tcrossprod(beta-t(betamean),beta-t(betamean))
    tausq=rgamma(1,shape=f/2,rate=betarms/2)
    tausq=1/tausq
    zrms=zmu-tcrossprod(x,beta)
    zrms=crossprod(zrms,zrms)
    sigsq=rgamma(1,shape=n/2,rate=zrms/2)
    sigsq=1/sigsq
    if(i>1000){
      sigsqsample[i-1000]=sigsq
    }
    muzero=rnorm(1,mean(mu),sqrt(lambdasq))
    murms=crossprod(mu-muzero,mu-muzero)
    lambdasq=rgamma(1,shape=(d+m)/2,rate=(eta+murms)/2)
    lambdasq=1/lambdasq
  }
  return=list(betasample=betasample,betavarsample=betavarsample,musample=musample,muvarsample=muvarsample,sigsqsample=sigsqsample)
}