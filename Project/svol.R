library(stats)
source('~/GitHub/SDS383D-course-work/Project/functions.R')

x=-106
sigma2=0.001
mu=1
#y=0.1
#sigma2=0.1
#mu=1
h=seq(from=1,to=10,by=0.01)
s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
r=(s-1)*exp(mu+sigma2/2)+x^2/2
qmode=r/(1+s)
c=p(x,qmode,sigma2,mu)/dinvgamma(qmode,scale=r,shape=s)
plot(p(x,h,sigma2,mu)~h,type='l',col="black",ylab="")
lines(c*dinvgamma(h,shape=s,scale=r)~h,col="blue")
lines(1.1*c*dinvgamma(h,shape=s,scale=r)~h,col="red")
lines(1.2*c*dinvgamma(h,shape=s,scale=r)~h,col="green")
legend(x=6.8,y=0.25,legend=c("p","c=1","c=1.1","c=1.2"),fill=c("black","blue","red","green"))


sp <- read.csv("~/GitHub/SDS383D-course-work/Project/SP1980-1987.csv")
ix<-read.csv("~/GitHub/SDS383D-course-work/Project/IXIC2007-2010 w.csv")
plot(ix$Close,type="l")
ar1=ar.ols(x=sp$Close,order.max=1)
lines(ar1$ar[1,1,1]*(sp$S.P.500-ar1$x.mean)+ar1$x.mean,col="red")
y=ar1$resid
y=na.omit(y)
plot(y,type='l')

n=length(ix$Close)
logchange=log(ix$Close[2:n])-log(ix$Close[1:(n-1)])
plot(logchange,type='l')


y1=(sp$Close[2:n]-sp$Close[1:(n-1)])/sp$Close[1:(n-1)]
plot(y1,type='l')

sample=sampler(logchange,1,1,0,10,0,10,1000,0,0.0000001)
plot(sample$delta_sample,type='l',ylab="",main="delta")
hist(sample$delta_sample,xlab='delta',ylab="",main="delta")
plot(sample$alpha_sample,type='l',ylab="",main="alpha")
hist(sample$alpha_sample,xlab='alpha',ylab="",main="alpha")
plot(sample$sigma_nu2_sample,type='l',ylab="",main="sigma_nu2")
hist(sample$sigma_nu2_sample,ylab="",main="sigma_nu2")
plot(log(sample$h_sample[,1000]),type='l',ylab="log(h)",main="1000th iteration")
plot(log(sample$h_sample[100,]),type='l',ylab="log(h)",main="evolution of h")
sample$rejectionrate
