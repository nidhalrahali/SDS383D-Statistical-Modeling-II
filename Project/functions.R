library(MCMCpack)
library(truncnorm)

p=function(y,h,sigma2,mu){
  return=h^(-1.5)*exp(-y^2/(2*h)-(log(h)-mu)^2/(2*sigma2))
}

nexth=function(y,h,sigma2,mu){
  s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
  r=(s-1)*exp(mu+sigma2/2)+y^2/2
  qmode=r/(1+s)
  c=1.1*p(y,qmode,sigma2,mu)/dinvgamma(qmode,scale=r,shape=s)
  hn=rinvgamma(1,shape=s,scale=r)
  accept=p(y,hn,sigma2,mu)/(c*dinvgamma(hn,scale=r,shape=s))
  roll=runif(1)
  if(roll>accept)hn=h
  return=hn
}

p_delta=function(delta,sigma_nu2,alpha){
  return=sqrt(1-delta^2)*exp(-alpha^2/((1-delta)*sigma_nu2))
}

nextdelta=function(delta,deltamean,deltasd,sigma_nu2,alpha){
  x=seq(from=-1,to=1,by=0.01)
  y=sqrt(1-x^2)*exp(-alpha^2/((1-x)*sigma_nu2))
  meanx=sum(x*y)/sum(y)
  varx=sum(x^2*y)/sum(y)
  varx=varx-meanx^2
  #print(deltasd)
  deltan=rtruncnorm(1,a=-1,b=1,mean=(deltamean*varx+meanx*deltasd^2)/(varx+deltasd^2),sd=deltasd*sqrt(varx)/sqrt(varx+deltasd^2))
  return=deltan
}

computemu=function(ln_h,delta,alpha){
  n=length(ln_h)
  mu=rep(0,n)
  mu[1]=delta*ln_h[2]+alpha
  mu[n]=delta*ln_h[n-1]+alpha
  for(i in 2:(n-1))mu[i]=(delta*(ln_h[i+1]+ln_h[i-1])+(1-delta)*alpha)/(1+delta^2)
  return=mu
}



sampler=function(y,nu_0,s_0,delta_0,sigma_delta2,alpha_0,sigma_alpha2,t){
  n=length(y)
  sigma_nu2=rinvgamma(1,shape=(nu_0)/2,scale=s_0/2)
  delta=rtruncnorm(1,a=-1,b=1,mean=delta_0,sd=sqrt(sigma_delta2))
  alpha=rnorm(1,mean=alpha_0,sd=sqrt(sigma_alpha2))
  h=sqrt(abs(y))
  h_sample=matrix(nrow=n,ncol=t)
  alpha_sample=rep(0,t)
  delta_sample=rep(0,t)
  sigma_nu2_sample=rep(0,t)
  for(ite in 1:t){
    ln_h=log(h)
    s1=sum(ln_h)
    s2=sum(ln_h^2)
    s1p=s1-ln_h[1]-ln_h[n]
    s2p=s2-ln_h[1]^2-ln_h[n]^2
    s3=0
    for(i in 1:(n-1))s3=s3+ln_h[i]*ln_h[i+1]
    sigma_nu_scale=(s_0+(n+2*delta/(1-delta))*alpha^2+s2+delta^2*s2p-2*alpha*(s1-delta*s1p)-2*delta*s3)/2
    sigma_nu2=rinvgamma(1,shape=(nu_0+n)/2,scale=sigma_nu_scale)
    deltamean=(sigma_nu2*delta_0+sigma_delta2*(s3-alpha*s1p))/(sigma_nu2+sigma_delta2*s2p)
    deltasd=sqrt(sigma_nu2*sigma_delta2/(sigma_nu2+sigma_delta2*s2p))
    delta=nextdelta(delta,deltamean,deltasd,sigma_nu2,alpha)
    #sigma_delta2=rinvgamma(1,shape=1/2,scale=(delta-delta_0)^2/2)
    alphamean=(sigma_alpha2*(s1-delta*s1p)+sigma_nu2*alpha_0)/(sigma_nu2+(n+2*delta/(1-delta))*sigma_alpha2)
    alpha=rnorm(1,mean=alphamean,sd=sqrt(sigma_nu2*sigma_alpha2/(sigma_nu2+(n+2*delta/(1-delta))*sigma_alpha2)))
    #sigma_alpha2=rinvgamma(1,shape=1/2,scale=(alpha-alpha_0)^2/2)
    mu=computemu(ln_h,delta,alpha)
    newh=h
    print(delta)
    newh[1]=nexth(y[1],h[1],sigma_nu2,mu[1])
    newh[n]=nexth(y[1],h[1],sigma_nu2,mu[n])
    sigma2=sigma_nu2/(1+delta^2)
    for(i in 2:(n-1))newh[i]=nexth(y[i],h[i],sigma2,mu[i])
    h=newh
    alpha_sample[ite]=alpha
    delta_sample[ite]=delta
    sigma_nu2_sample[ite]=sigma_nu2
    h_sample[,ite]=h
  }
  return=list(h_sample=h_sample,alpha_sample=alpha_sample,delta_sample=delta_sample,sigma_nu2_sample=sigma_nu2_sample)
}
