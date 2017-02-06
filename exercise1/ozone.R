library(mlbench)
ozone = data(Ozone, package='mlbench')

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator
A=solve(t(x) %*% x) %*% t(x)
betahat = A %*% y

# compute betacov from RMS
dy=y-x %*% betahat
sig=t(dy)%*%dy/(203-10)
sigma=diag(sig[1][1],203,203)
betacov=A%*%sigma%*%t(A)

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
lm1 = lm(y~x-1)
summary(lm1)
betacovlm = vcov(lm1)

# compare the diagonal elements of results
(sqrt(diag(betacov))-sqrt(diag(betacovlm)))/sqrt(diag(betacovlm))
