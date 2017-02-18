# this script generate a sample from multi-variate normal distribution,
# then use bootstrapp to get a sample distribution of estimated parameters

library(mvtnorm)
d = 2#dimension of distribution
n = 1000 #number of samples

# set parameter and generate sample x
mu = rep(2,d)
sigma = matrix(c(1,1,1,2),nrow=2,ncol=2)
x = rmvnorm(n,mu,sigma)
plot(x)

# get maximum likelyhood estimation of mu and sigma
mumle = colSums(x)/n
sigmamle = cov(x)
m = 800 # number of bootstrap samples

# do bootstrap
muboot = matrix(0, ncol=d, nrow=m)
ldet=0
for(i in 1:m){
  myindices = sample(1:n, n, replace=TRUE)
  bootx = subset(t(x),select=myindices)
  bootx = t(bootx)
  # for each sample, compute maximum-likelhood estimation of mean and covariance
  # also accumulate log(det(covariance))
  muboot[i,] = colSums(bootx)/n
  sigmaboot = cov(bootx)
  ldet=ldet+log(det(sigmaboot))
}

# check the distribution of muboot is a multivariate normal distribution
cov(muboot)# should be close to sigma/n
colSums(muboot)/m # should be close to mu
plot(muboot)

# check sigmaboot follows Wishart distribution by checking the formula
# of log expectations (https://en.wikipedia.org/wiki/Wishart_distribution)
ldet/m-digamma((n-1)/2)-digamma((n-2)/2)-log(det(sigma))-d*log(2/(n-1))
