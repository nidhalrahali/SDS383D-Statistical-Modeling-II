library(stats)
source('~/GitHub/SDS383D-course-work/Project/functions.R')

y=20
sigma2=0.2
mu=7
h=seq(from=0.1,to=100,by=0.1)
s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
r=(s-1)*exp(mu+sigma2/2)+y^2/2
qmode=r/(1+s)
c=p(y,qmode,sigma2,mu)/dinvgamma(qmode,scale=r,shape=s)
plot(p(y,h,sigma2,mu)~h,type='l',col="black",ylim=c(0,0.25),ylab="")
lines(c*dinvgamma(h,shape=s,scale=r)~h,col="blue")
lines(1.1*c*dinvgamma(h,shape=s,scale=r)~h,col="red")
lines(1.2*c*dinvgamma(h,shape=s,scale=r)~h,col="green")
legend(x=6.8,y=0.25,legend=c("p","c=1","c=1.1","c=1.2"),fill=c("black","blue","red","green"))


sp <- read.csv("~/GitHub/SDS383D-course-work/Project/SP.csv")
sp=sp[c(1,2,3,4)]
sp=sp[seq(from=1,to=252),]
plot(sp$S.P.500,type="l")
ar1=ar.ols(x=sp$S.P.500,order.max=1)
lines(ar1$ar*sp$S.P.500+ar1$x.intercept~seq(from=2,to=253),col="red")
y=ar1$resid
y=na.omit(y)
plot(y,type='l')
log(sd(y))

sample=sampler(y,1,1,0.5,0.1,1.2,0.1,1000)

