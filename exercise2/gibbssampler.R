library(mvtnorm)

gibbssampler=function(y,x,h,k,d,t){
  lambda=diag(1,dim(y)[1])
  omega=1
  result=matrix(nrow=2,ncol=t)
  for(i in 1:t){
    kstar=t(x)%*%lambda%*%x+k
    m=solve(kstar)%*%t(x)%*%lambda%*%y
    beta=rmvnorm(1,m,solve(omega*kstar))
    beta=t(beta)
    result[1,i]=beta[1]
    result[2,i]=beta[2]
    eta=t(y)%*%lambda%*%y-t(m)%*%kstar%*%m
    omega=rgamma(1,shape=d/2,rate=eta/2)
    alpha=y-x%*%beta
    for(j in 1:dim(y)[1]){
      lambda[j,j]=rgamma(1,shape=(h+1)/2,rate=h/2+alpha[j]*alpha[j]*omega/2)
    }
  }
  return=result
}

