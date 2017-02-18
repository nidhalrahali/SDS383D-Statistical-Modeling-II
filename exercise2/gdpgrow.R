# this script fit a linear heavy tail model of the relation between
# gdp growth rate and expense on defense 

# read data
mydata <- read.csv("~/courses/statistical modeling/gdpgrowth.csv", header=FALSE)
growth = as.matrix(mydata[2:dim(mydata)[1],3])
growth = as.numeric(growth)
defense = as.matrix(mydata[2:dim(mydata)[1],8])
defense = as.numeric(defense)

# plot data point
datapoint = cbind(defense,growth)
plot(datapoint)

#add intercept
defense = cbind(1,defense)

# set covariance for prior distribution
k = diag(1,2)
kstar = t(defense) %*% defense+k

# fit model
beta = solve(kstar)%*%t(defense)%*%growth

# plot result
abline(a = beta[1],b = beta[2])
