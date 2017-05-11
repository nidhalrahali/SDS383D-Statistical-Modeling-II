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

computemu=function(ln_h,delta,alpha){
  n=length(ln_h)
  mu=rep(0,n)
  mu[1]=(ln_h[2]-alpha)/delta
  mu[n]=delta*ln_h[n-1]+alpha
  for(i in 2:(n-1))mu[i]=(delta*(ln_h[i+1]+ln_h[i-1])+(1-delta)*alpha)/(1+delta^2)
  return=mu
}



sampler=function(y,nu_0,s_0,delta_0,sigma_delta2,alpha_0,sigma_alpha2,t,burnin){
  n=length(y)
  delta=1
  alpha=0
  h0=sd(y)^2
  h=rep(h0,n)
  h_sample=matrix(nrow=n,ncol=t)
  alpha_sample=rep(0,t)
  delta_sample=rep(0,t)
  sigma_nu2_sample=rep(0,t)
  for(ite in 1:(t+burnin)){
    if(ite%%100==0)print(ite)
    ln_h=log(h)
    s1=sum(ln_h)
    s2=sum(ln_h^2)
    s3=sum(ln_h[1:(n-1)]*ln_h[2:n])
    sigma_nu_scale=(s_0+(n-1)*alpha^2+(1+delta^2)*s2-ln_h[1]^2-delta^2*ln_h[n]^2-2*alpha*((1-delta)*s1-ln_h[1]+delta*ln_h[n])-2*delta*s3)/2
    sigma_nu2=rinvgamma(1,shape=(nu_0+n-1)/2,scale=sigma_nu_scale)
    
    deltamean=(sigma_nu2*delta_0+sigma_delta2*(s3-alpha*(s1-ln_h[n])))/(sigma_nu2+sigma_delta2*(s2-ln_h[n]^2))
    deltasd=sqrt(sigma_nu2*sigma_delta2/(sigma_nu2+sigma_delta2*(s2-ln_h[n]^2)))
    delta=rtruncnorm(1,a=-1,b=1,mean=deltamean,sd=deltasd)
    
    alphamean=(sigma_alpha2*((1-delta)*s1-ln_h[1]+delta*ln_h[n])+sigma_nu2*alpha_0)/(sigma_nu2+(n-1)*sigma_alpha2)
    alphasd=sqrt(sigma_nu2*sigma_alpha2/(sigma_nu2+(n-1)*sigma_alpha2))
    alpha=rnorm(1,mean=alphamean,sd=alphasd)
    
    mu=computemu(ln_h,delta,alpha)
    newh=h
    #print(delta)
    #newh[1]=nexth(y[1],h[1],sigma_nu2/delta^2,mu[1])
    #newh[n]=nexth(y[n],h[n],sigma_nu2,mu[n])
    sigma2=sigma_nu2/(1+delta^2)
    for(i in 2:(n-1))newh[i]=nexth(y[i],h[i],sigma2,mu[i])
    h=newh
    eps=rnorm(1)
    h[1]=exp(alpha+delta*log(h[2])+eps*sqrt(sigma_nu2))
    eps=rnorm(1)
    h[n]=exp(alpha+delta*log(h[n-1])+eps*sqrt(sigma_nu2))
    if(ite>burnin){
    alpha_sample[ite-burnin]=alpha
    delta_sample[ite-burnin]=delta
    sigma_nu2_sample[ite-burnin]=sigma_nu2
    h_sample[,ite-burnin]=h
  }
  }
  return=list(h_sample=h_sample,alpha_sample=alpha_sample,delta_sample=delta_sample,sigma_nu2_sample=sigma_nu2_sample)
}
