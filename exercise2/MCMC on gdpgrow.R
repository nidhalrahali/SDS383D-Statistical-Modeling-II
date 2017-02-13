# read data
mydata <- read.csv("~/courses/statistical modeling/gdpgrowth.csv", header=FALSE)
growth=as.matrix(mydata[2:dim(mydata)[1],3])
growth=as.numeric(growth)
growth=as.matrix(growth)
defense=as.matrix(mydata[2:dim(mydata)[1],8])
defense=as.numeric(defense)

# plot data point
datapoint=cbind(defense,growth)
plot(datapoint)

# add intercept
defense = cbind(1,defense)

k=diag(1,2)
d=80
beta=gibbssampler(growth,defense,4,k,d,5000)
plot(t(beta))
