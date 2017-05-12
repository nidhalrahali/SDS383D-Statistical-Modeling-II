library(MCMCpack)
library(truncnorm)

ln_p=function(y,h,sigma2,mu){
  return=-1.5*log(h)-y^2/(2*h)-(log(h)-mu)^2/(2*sigma2)
}

updateh=function(y,h,lhleft,lhright,alpha,delta,sigma_nu2){
  sigma2=sigma_nu2/(1+delta^2)
  mu=(delta*(lhleft+lhright)+(1-delta)*alpha)/(1+delta^2)
  s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
  r=(s-1)*exp(mu+sigma2/2)+y^2/2
  qmode=r/(1+s)
  ln_c=log(1.2)+ln_p(y,qmode,sigma2,mu)-log(dinvgamma(qmode,scale=r,shape=s))
  hn=rinvgamma(1,shape=s,scale=r)
  accept=exp(ln_p(y,hn,sigma2,mu)-ln_c-log(dinvgamma(hn,scale=r,shape=s)))
  roll=runif(1)
  rej=0
  if(roll>accept){
    hn=h
    rej=1
  }
  return=list(newh=hn,reject=rej)
}

p_ratio=function(y,h1,h2,sigma2,mu){
  return=(h2/h1)^1.5*exp(-y^2/(2*h1)+y^2/(2*h2)-((log(h1)-mu)^2-(log(h2)-mu)^2)/(2*sigma2))
}

updateh_rw=function(y,h,lhleft,lhright,alpha,delta,sigma_nu2,sd_rw){
  sigma2=sigma_nu2/(1+delta^2)
  mu=(delta*(lhleft+lhright)+(1-delta)*alpha)/(1+delta^2)
  hn=rtruncnorm(1,a=0,mean=h,sd=sd_rw)
  accept=p_ratio(y,hn,h,sigma2,mu)
  roll=runif(1)
  rej=0
  if(roll>accept){
    hn=h
    rej=1
  }
  return=list(newh=hn,reject=rej)
}

sampler=function(y,nu_0,s_0,delta_0,sigma_delta2,alpha_0,sigma_alpha2,t,burnin,sd_rw){
  n=length(y)
  delta=1
  alpha=0
  #h0=sd(y)^2
  #h=rep(h0,n)
  h=rep(10,n)
  ln_h=log(h)
  h_sample=matrix(nrow=n,ncol=t)
  alpha_sample=rep(0,t)
  delta_sample=rep(0,t)
  sigma_nu2_sample=rep(0,t)
  upd=0
  rej=0
  for(ite in 1:(t+burnin)){
    if(ite%%100==0)print(ite)
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
    #print(delta)
    
    for(i in 2:(n-1)){
        #r=updateh(y[i],h[i],ln_h[i-1],ln_h[i+1],delta,alpha,sigma_nu2)
        r=updateh_rw(y[i],h[i],ln_h[i-1],ln_h[i+1],delta,alpha,sigma_nu2,sd_rw)
        h[i]=r$newh
        rej=rej+r$reject
        upd=upd+1
        ln_h[i]=log(h[i])
    }
    ln_h[1]=rnorm(1,mean=alpha+delta*ln_h[2],sd=sqrt(sigma_nu2))
    h[1]=exp(ln_h[1])
    ln_h[n]=rnorm(1,mean=alpha+delta*ln_h[n-1],sd=sqrt(sigma_nu2))
    h[n]=exp(ln_h[n])
    
    if(ite>burnin){
    alpha_sample[ite-burnin]=alpha
    delta_sample[ite-burnin]=delta
    sigma_nu2_sample[ite-burnin]=sigma_nu2
    h_sample[,ite-burnin]=h
  }
  }
  return=list(h_sample=h_sample,alpha_sample=alpha_sample,delta_sample=delta_sample,sigma_nu2_sample=sigma_nu2_sample,rejectionrate=rej/upd)
}
