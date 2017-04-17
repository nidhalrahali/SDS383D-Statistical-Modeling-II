library(lattice)
library(dplyr)
droslong <- read.csv("~/GitHub/SDS383D-course-work/exercise4/gene/droslong.csv", header=TRUE)
droslong=droslong[c(1,2,4,5,6)]
droslong$group=as.numeric(droslong$group)
xyplot(log2exp~time|gene,data=droslong)
xyplot(log2exp~time|group,data=droslong%>%filter(replicate=="A"))
xyplot(log2exp~time|replicate,data=(droslong%>%filter(gene=="142798_at")))
