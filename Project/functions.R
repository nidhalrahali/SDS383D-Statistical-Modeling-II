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

computemu=function(ln_h,ln_h0,ln_hlast,delta,alpha){
  n=length(ln_h)
  mu=rep(0,n)
  mu[1]=(delta*(ln_h0+ln_h[2])+(1-delta)*alpha)/(1+delta^2)
  mu[n]=(delta*(ln_hlast+ln_h[n-1])+(1-delta)*alpha)/(1+delta^2)
  for(i in 2:(n-1))mu[i]=(delta*(ln_h[i+1]+ln_h[i-1])+(1-delta)*alpha)/(1+delta^2)
  return=mu
}

inith=function(y,ln_h0,ln_hlast,delta,alpha,sigma_nu2){
  n=length(y)
  ln_h=rep(0,n)
  for(i in 1:n)ln_h[i]=rnorm(1,mean=alpha/(1-delta),sd=sqrt(sigma_nu2/(1-delta^2)))
  mu=computemu(ln_h,ln_h0,ln_hlast,delta,alpha)
  h=exp(ln_h)
  newh=rep(0,n)
  sigma2=sigma_nu2/(1+delta^2)
  for(i in 1:n){
    while(1==1){
      newh[i]=nexth(y[i],h[i],sigma2,mu[i])
      if(newh[i]!=h[i])break
    }
  }
  return=newh
}
sampler=function(y,nu_0,s_0,delta_0,sigma_delta2,alpha_0,sigma_alpha2,t){
  n=length(y)
  sigma_nu2=rinvgamma(1,shape=(nu_0)/2,scale=s_0/2)
  delta=rnorm(1,mean=delta_0,sd=sqrt(sigma_delta2))
  alpha=rnorm(1,mean=alpha_0,sd=sqrt(sigma_alpha2))
  ln_h0=rnorm(1,mean=alpha/(1-delta),sd=sqrt(sigma_nu2/(1-delta^2)))
  ln_hlast=rnorm(1,mean=alpha/(1-delta),sd=sqrt(sigma_nu2/(1-delta^2)))
  h=inith(y,ln_h0,ln_hlast,delta,alpha,sigma_nu2)
  h_sample=matrix(nrow=n,ncol=t)
  alpha_sample=rep(0,t)
  delta_sample=rep(0,t)
  sigma_nu2_sample=rep(0,t)
  for(ite in 1:t){
    ln_h=log(h)
    sum_ln_h=sum(ln_h)
    sum_ln_h2=sum(ln_h^2)
    ln_product=ln_h0*ln_h[1]+ln_h[n]*ln_hlast
    for(i in 1:(n-1))ln_product=ln_product+ln_h[i]*ln_h[i+1]
    sigma_nu_scale=(s_0+(n+1)*alpha^2+(1+delta^2)*sum_ln_h2+delta^2*ln_h0^2+ln_hlast^2+2*alpha*(delta*ln_h0-(1-delta)*sum_ln_h-ln_hlast)-2*delta*ln_product)/2
    sigma_nu2=rinvgamma(1,shape=(nu_0+n+1)/2,scale=sigma_nu_scale)
    deltamean=(sigma_nu2*delta_0+sigma_delta2*(ln_product-alpha*(sum_ln_h+ln_h0)))/(sigma_nu2+sigma_delta2*(sum_ln_h2+ln_h0^2))
    delta=rtruncnorm(1,b=1,mean=deltamean,sd=sqrt(sigma_nu2*sigma_delta2/(sigma_nu2+sigma_delta2*(sum_ln_h2+ln_h0^2))))
    sigma_delta2=rinvgamma(1,shape=1/2,scale=(delta-delta_0)^2/2)
    alphamean=(sigma_alpha2*((1-delta)*sum_ln_h+ln_hlast-delta*ln_h0)+sigma_nu2*alpha_0)/(sigma_nu2+n*sigma_alpha2)
    alpha=rnorm(1,mean=alphamean,sd=sqrt(sigma_nu2*sigma_alpha2/(sigma_nu2+sigma_alpha2*n)))
    sigma_alpha2=rinvgamma(1,shape=1/2,scale=(alpha-alpha_0)^2/2)
    mu=computemu(ln_h,ln_h0,ln_hlast,delta,alpha)
    newh=h
    sigma2=sigma_nu2/(1+delta^2)
    print(ln_h[1])
    for(i in 1:n)newh[i]=nexth(y[i],h[i],sigma2,mu[i])
    h=newh
    eps=rnorm(1,mean=0,sd=1)
    ln_h0=(ln_h[1]-sqrt(sigma_nu2)*eps-alpha)/delta
    eps=rnorm(1,mean=0,sd=1)
    ln_hlast=alpha+delta*ln_h[n]+eps*sqrt(sigma_nu2)
    alpha_sample[ite]=alpha
    delta_sample[ite]=delta
    sigma_nu2_sample[ite]=sigma_nu2
    h_sample[,ite]=h
  }
  return=list(h_sample=h_sample,alpha_sample=alpha_sample,delta_sample=delta_sample,sigma_nu2_sample)
}
