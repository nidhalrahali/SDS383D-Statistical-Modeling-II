library(lattice)
library(dplyr)
source('~/GitHub/SDS383D-course-work/exercise4/gene/functions.R')
droslong <- read.csv("~/GitHub/SDS383D-course-work/exercise4/gene/droslong.csv", header=TRUE)
droslong=droslong[c(1,2,4,5,6)]
droslong$group=as.numeric(droslong$group)
droslong$genelabel=as.numeric(droslong$gene)
xyplot(log2exp~time|gene,data=droslong)
xyplot(log2exp~time|group,data=droslong%>%filter(replicate=="A"))
xyplot(log2exp~time|replicate,data=(droslong%>%filter(genelabel==4)))

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
sigma_eps_sample=sample$se
hist(sigma_eps_sample)

sigma_tau_sample=sample$st

par(mfrow=c(4,4))
for(i in 1:14)hist(sigma_tau_sample[i,],main=genename[i])

sigma_f_sample=sample$sf
par(mfrow=c(2,2))
for(i in 1:3)hist(sigma_f_sample[i,],main=i)

sigma_h_sample=sample$sh
par(mfrow=c(4,4))
for(i in 1:14)hist(sigma_h_sample[i,],main=genename[i])

par(mfrow=c(2,2))
for(i in 1:3){
  groupsample=matrix(nrow=12,ncol=t)
  for(j in 1:t)groupsample[,j]=fsample[j,,i]
  boxplot(t(groupsample),xlab="time",ylab="log2exp",ylim=c(5,18),main=i)
}

genegroup=rep(0,14)
genename=rep(0,14)
for(i in 1:14){
  genegroup[gene[i]]=gr[i]
  genename[gene[i]]=as.character(droslong$gene[i])
}

par(mfrow=c(4,4))
for(i in 1:14){
  genesample=matrix(nrow=12,ncol=t)
  for(j in 1:t)genesample[,i]=hsample[j,,i]
  boxplot(t(genesample),xlab="time",ylab="log2exp",main=genename[i])
}

par(mfrow=c(3,5))
for(i in 1:14){
  genesample=matrix(nrow=12,ncol=t)
  for(j in 1:t)genesample[,i]=hsample[j,,i]+fsample[j,,genegroup[i]]
  boxplot(t(genesample),xlab="time",ylab="log2exp",ylim=c(5,18),main=genename[i])
  lines(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="A")),col="red")
  lines(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="B")),col="purple")
  lines(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="C")),col="blue")
}


