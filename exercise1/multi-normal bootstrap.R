library(mvtnorm)
d=2#dimension of distribution
n=1000 #number of samples
# set parameter and generate sample x
mu=rep(2,d)
sigma=matrix(c(1,1,1,2),nrow=2,ncol=2)
x=rmvnorm(n,mu,sigma)
plot(x)
# get mle of mu and sigma
mumle=colSums(x)/n
sigmamle=cov(x)
m=800 # number of bootstrap samples
# bootstrap
muboot = matrix(0, ncol=d, nrow=m)
ldet=0
for(i in 1:m){
  myindices = sample(1:n, n, replace=TRUE)
  bootx=subset(t(x),select=myindices)
  bootx=t(bootx)
  muboot[i,]=colSums(bootx)/n
  sigmaboot=cov(bootx)
  ldet=ldet+log(det(sigmaboot))
}
# check the distribution of muboot
cov(muboot)
colSums(muboot)/m

# check sigmaboot follows Wishart distribution
ldet/m-digamma((n-1)/2)-digamma((n-2)/2)-log(det(sigma))-d*log(2/(n-1))
