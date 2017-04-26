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
Cinv=solve(C)
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

genegroup=rep(0,14)
genename=rep(0,14)
for(i in 1:14){
  genegroup[gene[i]]=gr[i]
  genename[gene[i]]=as.character(droslong$gene[i])
}

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
  boxplot(t(groupsample),xlab="time",ylab="log2exp",ylim=c(5,20),main=i)
}


par(mfrow=c(4,4))
for(i in 1:14){
  genesample=matrix(nrow=12,ncol=t)
  for(j in 1:t)genesample[,j]=hsample[j,,i]
  boxplot(t(genesample),xlab="time",ylab="log2exp",,ylim=c(-8,8),main=genename[i])
}

par(mfrow=c(3,5))
for(i in 1:14){
  genesample=matrix(nrow=12,ncol=t)
  for(j in 1:t)genesample[,j]=hsample[j,,i]+fsample[j,,genegroup[i]]
  boxplot(t(genesample),xlab="time",ylab="log2exp",ylim=c(5,18),main=genename[i])
  lines(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="A")),col="red")
  lines(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="B")),col="purple")
  lines(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="C")),col="blue")
}

g=12
genesample=matrix(nrow=12,ncol=t)
for(j in 1:t)genesample[,j]=hsample[j,,g]+fsample[j,,genegroup[g]]
boxplot(t(genesample),xlab="time",ylab="log2exp",main=genename[g])
lines(log2exp~time,data=(droslong%>%filter(genelabel==g,replicate=="A")),col="red")
lines(log2exp~time,data=(droslong%>%filter(genelabel==g,replicate=="B")),col="purple")
lines(log2exp~time,data=(droslong%>%filter(genelabel==g,replicate=="C")),col="blue")

mygrid=seq(from=0.9,to=12.1,by=0.11)
Cstarstar=matrix(nrow=102,ncol=102)
for(i in 1:102){
  for(j in 1:102){
    Cstarstar[i,j]=materncov(abs(mygrid[i]-mygrid[j]),b)
  }
}
Cstar=matrix(nrow=102,ncol=12)
for(i in 1:102){
  for(j in 1:12){
    Cstar[i,j]=materncov(abs(mygrid[i]-j),b)
  }
}
Cdiag=diag(Cstarstar-Cstar%*%tcrossprod(Cinv,Cstar))
par(mfrow=c(3,5))
for(i in 1:14){
  hmean=matrix(0,nrow=12,ncol=1)
  fmean=matrix(0,nrow=12,ncol=1)
  for(j in 1:t){
    hmean=hmean+hsample[j,,i]
    fmean=fmean+fsample[j,,genegroup[i]]
  }
  hmean=hmean/t
  fmean=fmean/t
  pred=Cstar%*%Cinv%*%fmean+Cstar%*%Cinv%*%hmean
  sigma_f_inv_mean=mean(1/sigma_f_sample[genegroup[i],])
  sigma_h_inv_mean=mean(1/sigma_h_sample[i,])
  sigma_eps_inv_mean=mean(1/sigma_eps_sample)
  sigma_tau_inv_mean=mean(1/sigma_tau_sample[i,])
  prederror=1.96*sqrt((1/sigma_f_inv_mean+1/sigma_h_inv_mean)*Cdiag+1/sigma_eps_inv_mean+1/sigma_tau_inv_mean)
  plot(pred~mygrid,xlab="time",ylab="log2exp",ylim=c(5,18),type="l",main=genename[i])
  polygon(x=c(mygrid,rev(mygrid)),y=c(pred+prederror,rev(pred-prederror)),col="lightgrey")
  lines(pred~mygrid)
  points(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="A")),col="red")
  points(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="B")),col="purple")
  points(log2exp~time,data=(droslong%>%filter(genelabel==i,replicate=="C")),col="blue")
}