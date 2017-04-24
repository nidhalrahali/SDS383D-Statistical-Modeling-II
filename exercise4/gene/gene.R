library(lattice)
library(dplyr)
source('~/GitHub/SDS383D-course-work/exercise4/gene/functions.R')
droslong <- read.csv("~/GitHub/SDS383D-course-work/exercise4/gene/droslong.csv", header=TRUE)
droslong=droslong[c(1,2,4,5,6)]
droslong$group=as.numeric(droslong$group)
droslong$genelabel=as.numeric(droslong$gene)
xyplot(log2exp~time|gene,data=droslong)
xyplot(log2exp~time|group,data=droslong%>%filter(replicate=="A"))
xyplot(log2exp~time|replicate,data=(droslong%>%filter(genelabel==1)))

b=1
C=matrix(nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    C[i,j]=materncov(abs(i-j),b)
  }
}
gene=droslong$genelabel
y=droslong$log2exp
gr=droslong$group
ti=droslong$time
t=1000

sample=gibbssampler(y,gene,gr,ti,C,t)
hsample=sample$h
fsample=sample$f
genesample=matrix(nrow=12,ncol=t)
for(i in 1:t)genesample[,i]=fsample[i,,1]+hsample[i,,3]
boxplot(t(genesample),xlab="time",ylab="log2exp")
lines(log2exp~time,data=(droslong%>%filter(genelabel==1,replicate=="A")))
lines(log2exp~time,data=(droslong%>%filter(genelabel==1,replicate=="B")))
lines(log2exp~time,data=(droslong%>%filter(genelabel==1,replicate=="C")))
