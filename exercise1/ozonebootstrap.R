# this script estimate the covariance matrix of the estimator (for the
# linear model of ozone concentration) using bootstrap. 
library(mlbench)
ozone = data(Ozone, package='mlbench')

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = as.matrix(ozone[,1])
x = as.matrix(ozone[,2:10])
# add an intercept
x = cbind(1,x)
# choose the number of samples in bootstrap
n=500
# bootstrap
boot = matrix(0, ncol=10, nrow=n)
for(i in 1:n){
  myindices = sample(1:nrow(y), nrow(y), replace=TRUE)
  lmtemp = lm(y~x-1,subset=myindices)
  boot[i,] = coef(lmtemp)
}
betacovboot=cov(boot)

# compute the cov using parametric estimation
A=solve(t(x) %*% x) %*% t(x)
betahat = A %*% y
dy=y-x %*% betahat
sigma=t(dy)%*%dy/(203-10)
betacov=A%*%t(A)*sigma[1][1]

# compute cov using lm function
lm1 = lm(y~x-1)
betacovlm = vcov(lm1)

# compare the diagonal elements of covariance matrix
(diag(betacovboot)-diag(betacovlm))/diag(betacovlm)
(diag(betacovboot)-diag(betacov))/diag(betacov)

#compare all elements of covariance matrix, the off diagonal elements do not match well
(betacovboot-betacovlm)/betacovlm
(betacovboot-betacov)/betacovlm
