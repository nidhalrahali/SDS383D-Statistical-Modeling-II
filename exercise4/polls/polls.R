library(lme4)
library(dplyr)#$filter
library(mosaic)
polls <- read.csv("~/GitHub/SDS383D-course-work/exercise4/polls/polls.csv")
polls=na.omit(polls)
polls=polls[c(4,5,6,7,8,9)]
summary(polls)

edulevel=c('NoHS','HS','SomeColl','Bacc')
pollbyedu=rep(0,length(edulevel))
educount=pollbyedu
for (i in 1:length(edulevel)){
  educount[i]=count(~edu==edulevel[i],data=polls)
  pollbyedu[i]=perc(~bush==1,data=polls%>%filter(edu==edulevel[i]))
}
plot(pollbyedu,ylab="support",type='h')
plot(educount,ylab="number of polls",type='h')

agelevel=c("18to29","30to44","45to64","65plus")
pollbyage=rep(0,length(agelevel))
agecount=pollbyage
for (i in 1:length(agelevel)){
  agecount[i]=count(~age==agelevel[i],data=polls)
  pollbyage[i]=perc(~bush==1,data=polls%>%filter(age==agelevel[i]))
}
plot(agecount,ylab="number of polls",type='h')
plot(pollbyage,ylab="upport",type='h')

statelist=levels(polls$state)
statecount=rep(0,length(statelist))
pollbystate=statecount
for(i in 1:length(statelist)){
  statecount[i]=count(~state==statelist[i],data=polls)
  pollbystate[i]=perc(~bush==1,data=polls%>%filter(state==statelist[i]))
}
plot(pollbystate,ylab="bush support rate",type='h')
plot(statecount,ylab="bush support rate",type='h')
plot(pollbystate~statecount,xlab="number of polls",ylab="bush support rate")

polls$edu=factor(polls$edu,levels=edulevel)
polls$age=factor(polls$age,levels=agelevel)

model=lmer(bush~edu+age+female+black+(1|state),data=polls)
summary(model)

ranef(model)
plot(ranef(model)$state[,1]~pollbystate)


polls$pred=predict(model,polls)
plot(polls$pred)
hist(polls$pred)
edupredict=rep(0,length(edulevel))
for (i in 1:length(edulevel)){
  edupredict[i]=perc(~pred>0.5,data=polls%>%filter(edu==edulevel[i]))
}
plot(edupredict)

agepredict=rep(0,length(agelevel))
for (i in 1:length(edulevel)){
  agepredict[i]=perc(~pred>0.5,data=polls%>%filter(age==agelevel[i]))
}
plot(agepredict)

statepredict=rep(0,length(statelist))
for(i in 1:length(statelist)){
  statepredict[i]=perc(~pred>0.5,data=polls%>%filter(state==statelist[i]))
}
plot(statepredict)
plot(statepredict~pollbystate)

plot(polls$pred~polls$weight)

state=as.numeric(polls$state)
y=as.vector(polls$bush)
x=matrix(0,nrow=length(y),ncol=9)
for(i in 1:length(y)){
  x[i,1]=1
  x[i,8]=polls$female[i]
  x[i,9]=polls$black[i]
  if(polls$edu[i]=="HS"){
    x[i,2]=1
  }
  if(polls$edu[i]=="SomeColl"){
    x[i,3]=1
  }
  if(polls$edu[i]=="Bacc"){
    x[i,4]=1
  }
  if(polls$age[i]=="65plus"){
    x[i,7]=1
  }
  if(polls$age[i]=="30to44"){
    x[i,5]=1
  }
  if(polls$age[i]=="45to64"){
    x[i,6]=1
  }
}
beta=rep(0,9)
mu=rep(0,49)
lambdasq=1
sigsq=1
d=1
eta=1
t=2000
tausq=1
gibbssampler(y,x,state,beta,mu,statecount,lambdasq,sigsq,tausq,d,eta,t)
