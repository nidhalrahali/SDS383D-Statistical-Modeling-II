library(MCMCpack)

p=function(y,h,sigma2,mu){
  return=h^(-1.5)*exp(-y^2/(2*h)-(log(h)-mu)^2/(2*sigma2))
}

nexth=function(y,h,sigma2,mu){
  s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
  r=(s-1)*exp(mu+sigma2/2)+y^2/2
  qmode=r/(1+s)
  c=1.1*p(y,qmode,sigma2,mu)/dinvgamma(qmode,scale=r,shape=s)
  hn=rinvgamma(1,shape=s,scale=r)
  accept=p(y,hn,sigma2,mu)/dinvgamma(hn,scale=r,shape=s)
  roll=runif(1)
  if(roll>accept)hn=h
  return=hn
}

sampler=function(y){
  for(ite in 1:t){
    
  }
}
