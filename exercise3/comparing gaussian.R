x=seq(0,1,by=0.001)

# same d,t2, different t1
g=as.vector(gaussian(materncov,x,0.05,10,0.001))
plot(g~x,type="l")

g=as.vector(gaussian(materncov,x,0.05,5,0.001))
lines(x,g,col="pink")

g=as.vector(gaussian(materncov,x,0.05,1,0.001))
lines(x,g,col="red")

g=as.vector(gaussian(materncov,x,0.05,0.1,0.001))
lines(x,g,col="blue")

# same t1,t2, different d
g=as.vector(gaussian(materncov,x,0.05,1,0.001))
plot(g~x,type="l")

g=as.vector(gaussian(materncov,x,0.5,1,0.001))
lines(x,g,col="pink")

g=as.vector(gaussian(materncov,x,1,1,0.001))
lines(x,g,col="red")

g=as.vector(gaussian(materncov,x,5,1,0.001))
lines(x,g,col="blue")

#t2=0
g=as.vector(gaussian(materncov,x,0.05,10,0))
plot(g~x,type="l")

# same parameters, different cov fuction
g=as.vector(gaussian(materncov,x,0.1,10,0.001))
plot(g~x,type="l")

g=as.vector(gaussian(secov,x,0.1,10,0.001))
lines(x,g,col="red")
