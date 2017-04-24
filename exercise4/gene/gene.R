library(lattice)
library(dplyr)
source('~/GitHub/SDS383D-course-work/exercise4/gene/functions.R')
droslong <- read.csv("~/GitHub/SDS383D-course-work/exercise4/gene/droslong.csv", header=TRUE)
droslong=droslong[c(1,2,4,5,6)]
droslong$group=as.numeric(droslong$group)
droslong$genelabel=as.numeric(droslong$gene)
xyplot(log2exp~time|gene,data=droslong)
xyplot(log2exp~time|group,data=droslong%>%filter(replicate=="A"))
xyplot(log2exp~time|replicate,data=(droslong%>%filter(gene=="142798_at")))

b=1
C=matrix(nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    C[i,j]=materncov(abs(i-j),b)
  }
}
Cinv=solve(10*C)
gene=droslong$genelabel
y=droslong$log2exp
gr=droslong$group
ti=droslong$time
