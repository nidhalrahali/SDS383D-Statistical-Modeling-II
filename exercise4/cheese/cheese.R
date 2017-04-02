library(lme4)
library(dplyr)
library(ggplot2)
library(mosaic)
cheese <- read.csv("~/GitHub/SDS383D-course-work/exercise4/cheese/cheese.csv")
summary(cheese)

cheese$store=as.integer(cheese$store)
cheese$store=factor(cheese$store)
cheese$logprice=log(cheese$price)
cheese$logvol=log(cheese$vol)


disp0=cheese %>% filter(disp == 0)
disp1=cheese %>% filter(disp==1)
disp0.mu = mean(disp0$logvol)
disp1.mu = mean(disp1$logvol)
p=ggplot(cheese, aes(disp,logvol)) + geom_jitter() + geom_hline(yintercept = mean(cheese$logvol),col="red")
p=p +  geom_hline(yintercept = disp0.mu ,col="green") + annotate("text", x = 1, y=12, label="Disp = 0, 1968 obs", col="green")
p=p +  geom_hline(yintercept = disp1.mu ,col="blue") + annotate("text", x = 2, y=12, label="Disp = 1, 3587 obs", col="blue")
p

xyplot(logvol ~ logprice | store, data=cheese, type = c("p", "r"),  group = disp, auto.key = list(), par.strip.text=list(cex=0.5))

plot(cheese$logvol~cheese$logprice)
boxplot(cheese$logvol~cheese$store)
boxplot(cheese$logvol~cheese$disp)
boxplot(cheese$logprice~cheese$disp)

model1=lmer(logvol~(disp+logprice+disp:logprice|store),data=cheese)
summary(model1)
ranef(model1)
plot(model1)

cheese$disp=factor(cheese$disp)
model2=lmer(logvol~logprice+(1+logprice|store:disp),data=cheese)
plot(coef(model2))
summary(model2)
ranef(model2)
plot(model2)
