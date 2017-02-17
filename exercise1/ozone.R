# this script builds a liner model of data of ozone concentration in LA on other atmospherical variable. We estimate
# the covariance matrix of estimator beta through RMS of residual
library(mlbench)
ozone = data(Ozone, package='mlbench')

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator betahat
A = solve(t(x) %*% x) %*% t(x)
betahat = A %*% y

# compute betacov from RMS
dy = y-x %*% betahat

# sig is the variance of error, from which we will compute the covariance of estimator
sigma = t(dy)%*%dy/(203-10)
betacov = A%*%t(A)*sigma[1][1]

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
lm1 = lm(y~x-1)
summary(lm1)
betacovlm = vcov(lm1)

# compare the diagonal elements of results by computing relative difference
(sqrt(diag(betacov))-sqrt(diag(betacovlm)))/sqrt(diag(betacovlm))
