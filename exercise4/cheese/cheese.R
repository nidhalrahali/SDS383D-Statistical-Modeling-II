library(lme4)
library(dplyr)#$filter
library(ggplot2)#ggplot
library(mosaic)#xyplot
source('~/GitHub/SDS383D-course-work/exercise4/cheese/gibbssampler.R')
cheese <- read.csv("~/GitHub/SDS383D-course-work/exercise4/cheese/cheese.csv")
summary(cheese)

cheese$store=as.integer(cheese$store)
cheese$store=factor(cheese$store)
cheese$logprice=log(cheese$price)
cheese$logvol=log(cheese$vol)


plot(cheese$logvol~cheese$logprice)
boxplot(cheese$logvol~cheese$store)
boxplot(cheese$logvol~cheese$disp)
boxplot(cheese$logprice~cheese$disp)

#visualize the general effect that adds has on the sale
disp0=cheese %>% filter(disp == 0)
disp1=cheese %>% filter(disp==1)
disp0.mu = mean(disp0$logvol)
disp1.mu = mean(disp1$logvol)
p=ggplot(cheese, aes(disp,logvol)) + geom_jitter() + geom_hline(yintercept = mean(cheese$logvol),col="red")
p=p +  geom_hline(yintercept = disp0.mu ,col="green") + annotate("text", x = 1, y=12, label="Disp = 0, 1968 obs", col="green")
p=p +  geom_hline(yintercept = disp1.mu ,col="blue") + annotate("text", x = 2, y=12, label="Disp = 1, 3587 obs", col="blue")
p

xyplot(logvol ~ logprice | store, data=cheese, type = c("p", "r"),  group = disp, auto.key = list(), par.strip.text=list(cex=0.5))


model1=lmer(logvol~(disp+logprice+disp:logprice|store),data=cheese)
summary(model1)
coef(model1)
ranef(model1)
boxplot(resid(model1) ~ store, data=cheese,  main='residuals by store') 
cheese$p1=predict(model1,cheese)
xyplot(p1 ~ logprice | store, data=cheese, type = c("p", "r"),  group = disp, auto.key = list(), par.strip.text=list(cex=0.5))

cheese$disp=factor(cheese$disp)
model2=lmer(logvol~(1+logprice|store:disp),data=cheese)
summary(model2)
coef(model2)
ranef(model2)
boxplot(resid(model2) ~ store:disp, data=cheese,  main='residuals by store:disp') 
cheese$p2=predict(model2,cheese)
#seems that it splits data too much and gives some weird result in the groups that lacks data
xyplot(p2 ~ logprice | store, data=cheese, type = c("p", "r"),  group = disp, auto.key = list(), par.strip.text=list(cex=0.5))

model3=lmer(logvol~(1+logprice|store)+(1+logprice|disp),data=cheese)
summary(model3)
coef(model3)
cheese$p3=predict(model3,cheese)
#almost no distinction between disp and non-disp 
xyplot(p3 ~ logprice | store, data=cheese, type = c("p", "r"),  group = disp, auto.key = list(), par.strip.text=list(cex=0.5))

cheese$pd=cheese$logprice*cheese$disp
xx=array(dim=c(88,4,4))
xy=matrix(nrow=4,ncol=88)
yy=rep(1,88)
for(i in 1:88){
  storedata=cheese%>%filter(store==i)
y=as.matrix(storedata$logvol)
x=as.matrix(cbind(1,storedata$logprice,storedata$disp,storedata$pd))
  xx[i,,]=crossprod(x,x)
  xy[,i]=crossprod(x,y)
  yy[i]=crossprod(y,y)
}
t=1000
sampling=gibbssampler(xx,xy,yy,t,5555)
betasample=sampling$betasample
sigmasample=sampling$sigmasample
siginvsample=sampling$siginvsample

hist(1/siginvsample,xlab="sigma^2")

betameanbystore=matrix(nrow=88,ncol=4)
for(i in 1:88){
  betameanbystore[i,]=colMeans(betasample[i,,])
}

lme4coef=as.matrix(coef(model1)$store)

plotstore=function(st){
  storedata=cheese%>%filter(store==st)
  xl=c(min(storedata$logprice),max(storedata$logprice))
  yl=c(min(storedata$logvol),max(storedata$logvol))
plot(logvol~logprice,data=storedata%>%filter(disp==0),xlim=xl,ylim=yl,col='blue')
points(logvol~logprice,data=storedata%>%filter(disp==1),col='red')
abline(betameanbystore[st,1],betameanbystore[st,2],col='blue')
abline(betameanbystore[st,1]+betameanbystore[st,3],betameanbystore[st,2]+betameanbystore[st,4],col='red')
abline(lme4coef[st,4],lme4coef[st,2],lty='dashed',col='blue')
abline(lme4coef[st,4]+lme4coef[st,1],lme4coef[st,3]+lme4coef[st,2],lty='dashed',col='red')
}
plotstore(15)
plotstore(30)
plotstore(40)

intercept=matrix(nrow=88,ncol=t)
slope=matrix(nrow=88,ncol=t)
dispeffect1=matrix(nrow=88,ncol=t)
dispeffect2=matrix(nrow=88,ncol=t)
for(i in 1:88){
  intercept[i,]=betasample[i,,1]
  slope[i,]=betasample[i,,2]
  dispeffect1[i,]=betasample[i,,3]
  dispeffect2[i,]=betasample[i,,4]
}
boxplot(t(dispeffect1))
points(lme4coef[,1],col="red")

boxplot(t(dispeffect2))
points(lme4coef[,3],col='red')

boxplot(t(intercept))
boxplot(t(slope))

cheese$betameanpredict=cheese$logvol
for(i in 1:5555){
  cheese$betameanpredict[i]=betameanbystore[cheese$store[i],]%*%c(1,cheese$logprice[i],cheese$disp[i],cheese$pd[i])
}
xyplot(betameanpredict ~ logprice | store, data=cheese, type = c("p"),  group = disp, auto.key = list(), par.strip.text=list(cex=0.5))
