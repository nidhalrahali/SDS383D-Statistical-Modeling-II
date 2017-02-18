# this script fit the linear model on the relation between the gdp
# growth rate and defense expanse by generating estimator using markov chain monte carlo 
# from a posterrior distribution
# read data
mydata <- read.csv("~/courses/statistical modeling/gdpgrowth.csv", header=FALSE)
growth = as.matrix(mydata[2:dim(mydata)[1],3])
growth = as.numeric(growth)
growth = as.matrix(growth)
defense = as.matrix(mydata[2:dim(mydata)[1],8])
defense = as.numeric(defense)

# plot data point
datapoint = cbind(defense,growth)
plot(datapoint)

# add intercept
defense = cbind(1,defense)

# set prior distribution parameters
k = diag(1,2)
d = 80

# source the function which does the sampling
source('~/GitHub/SDS383D-course-work/exercise2/gibbssampler.R')
beta = gibbssampler(growth,defense,4,k,d,5000)
# plot the result
abline(mean(beta[2]),mean(beta[1]))
