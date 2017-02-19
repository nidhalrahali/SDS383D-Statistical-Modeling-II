# gibbs sampler
library(mvtnorm)

# y,x: the vector of observables
# h, two times the shape of prior gamma distribution of lambda
# k: the prior covariance matrix
# d: two times the shape of prior gamma distribution of omega
# t: number of samples to be generated
# the return of this function is all the beta generated in the sample
# for further explanations of variables, see the excercise 2 solution.pdf
gibbssampler=function(y,x,h,k,d,t){
  lambda = diag(1,dim(y)[1])
  omega = 1
  result = matrix(nrow=2,ncol=t)
  for(i in 1:t+1000){
    kstar = t(x)%*%lambda%*%x+k
    m = solve(kstar)%*%t(x)%*%lambda%*%y
    beta = rmvnorm(1,m,solve(omega*kstar))
    beta = t(beta)
    if(i>1000){
    result[1,i-1000] = beta[1]
    result[2,i-1000] = beta[2]
    }
    eta = t(y)%*%lambda%*%y-t(m)%*%kstar%*%m
    omega = rgamma(1,shape=d/2,rate=eta/2)
    alpha = y-x%*%beta
    for(j in 1:dim(y)[1]){
      lambda[j,j] = rgamma(1,shape=(h+1)/2,rate=h/2+alpha[j]*alpha[j]*omega/2)
    }
  }
  return = result
}

