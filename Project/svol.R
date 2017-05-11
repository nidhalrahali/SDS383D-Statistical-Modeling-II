library(stats)
source('~/GitHub/SDS383D-course-work/Project/functions.R')

y=-0.00034
sigma2=0.0005
mu=-2.6
#y=0.1
#sigma2=0.1
#mu=1
h=seq(from=1,to=10,by=0.01)
s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
r=(s-1)*exp(mu+sigma2/2)+y^2/2
qmode=r/(1+s)
c=p(y,qmode,sigma2,mu)/dinvgamma(qmode,scale=r,shape=s)
plot(p(y,h,sigma2,mu)~h,type='l',col="black",ylab="")
lines(c*dinvgamma(h,shape=s,scale=r)~h,col="blue")
lines(1.1*c*dinvgamma(h,shape=s,scale=r)~h,col="red")
lines(1.2*c*dinvgamma(h,shape=s,scale=r)~h,col="green")
legend(x=6.8,y=0.25,legend=c("p","c=1","c=1.1","c=1.2"),fill=c("black","blue","red","green"))


sp <- read.csv("~/GitHub/SDS383D-course-work/Project/SP.csv")
plot(sp$S.P.500,type="l")
ar1=ar.ols(x=sp$S.P.500,order.max=1)
lines(ar1$ar[1,1,1]*(sp$S.P.500-ar1$x.mean)+ar1$x.mean,col="red")
y=ar1$resid
y=na.omit(y)
plot(y,type='l')

n=length(sp$S.P.500)
y1=(sp$S.P.500[2:n]-sp$S.P.500[1:(n-1)])/sp$S.P.500[1:(n-1)]
plot(y1,type='l')
ar2=ar.ols(x=y1,order.max=1)

sample=sampler(y1,1,1,0,10,0,10,2000,4000)
plot(sample$delta_sample,type='l',ylab="",main="delta")
plot(sample$alpha_sample,type='l',ylab="",main="alpha")
plot(sample$sigma_nu2_sample,type='l',ylab="",main="sigma_nu2")
plot(log(sample$h_sample[,1000]),type='l',ylab="log(h)",main="1000th iteration")
plot(log(sample$h_sample[100,]),type='l',ylab="log(h)",main="evolution of h")
