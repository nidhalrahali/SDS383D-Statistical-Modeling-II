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
p=rep(1,88)
for(i in 1:88){
  storedata=cheese%>%filter(store==i)
y=as.matrix(storedata$logvol)
p[i]=nrow(y)
x=as.matrix(cbind(1,storedata$logprice,storedata$disp,storedata$pd))
  xx[i,,]=crossprod(x,x)
  xy[,i]=crossprod(x,y)
  yy[i]=crossprod(y,y)
}
gibbssampler(xx,xy,yy,p,1000)
