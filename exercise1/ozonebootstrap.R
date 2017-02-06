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
sig=t(dy)%*%dy/(203-10)
sigma=diag(sig[1][1],203,203)
betacov=A%*%sigma%*%t(A)

# compute cov using lm function
lm1 = lm(y~x-1)
betacovlm = vcov(lm1)

# compare the diagonal elements of results
(diag(betacovboot)-diag(betacovlm))/diag(betacovlm)
(diag(betacovboot)-diag(betacov))/diag(betacov)

#compare all the results
(betacovboot-betacovlm)/betacovlm
(betacovboot-betacov)/betacovlm