library(lme4)
cheese <- read.csv("~/GitHub/SDS383D-course-work/exercise4/cheese/cheese.csv")
summary(cheese)

cheese$store=as.integer(cheese$store)
cheese$store=factor(cheese$store)
cheese$logprice=log(cheese$price)
cheese$logvol=log(cheese$vol)

plot(cheese$logvol~cheese$logprice)
boxplot(cheese$logvol~cheese$store)
boxplot(cheese$logvol~cheese$disp)

model1=lmer(logvol~logprice+disp+(1+disp+logprice|store),data=cheese)
summary(model1)
ranef(model1)
plot(model1)

cheese$disp=factor(cheese$disp)
model2=lmer(logvol~logprice+(1+logprice|store:disp),data=cheese)
summary(model2)
ranef(model2)
plot(model2)