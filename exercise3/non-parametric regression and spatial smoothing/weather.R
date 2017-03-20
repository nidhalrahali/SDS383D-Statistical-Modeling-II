source('~/GitHub/SDS383D-course-work/exercise3/gaussian process.R')
weather <- read.csv("~/GitHub/SDS383D-course-work/exercise3/weather.csv")
weather=as.matrix(weather)

temp=t(as.matrix(weather[,2]))
pres=t(as.matrix(weather[,1]))
lon=weather[,3]
lat=weather[,4]
n=length(lon)

b=30
t1=10
t2=0.001
sigsquare=0.01
b1=1
b2=1

C=matrix(ncol=n,nrow=n)
for(i in 1:n){
  for(j in i:n){
    d=sqrt((lon[i]-lon[j])^2+(lat[i]-lat[j])^2)
    C[i,j]=materncov(distance(lon[i],lat[i],lon[j],lat[j],b1,b2),b,t1,t2)
    C[j,i]=C[i,j]
  }
}
D=solve(C+sigsquare*diag(ncol(C)))
cstarstar=materncov(0,b,t1,t2)

longrid=t(as.matrix(seq(from=min(lon),to=max(lon),by=0.1)))
latgrid=t(as.matrix(seq(from=min(lat),to=max(lat),by=0.1)))
temphat=matrix(ncol=ncol(latgrid),nrow=ncol(longrid))
preshat=temphat
variance=temphat
for(i in 1:ncol(longrid)){
  for(j in 1:ncol(latgrid)){
  cstar=temp
  for(k in 1:n){
    cstar[k]=materncov(distance(longrid[i],latgrid[j],lon[k],lat[k],b1,b2),b,t1,t2)
  }
  W=tcrossprod(cstar,D)
  temphat[i,j]=tcrossprod(temp,W)
  preshat[i,j]=tcrossprod(pres,W)
  variance[i,j]=sqrt((cstarstar-tcrossprod(cstar,W)))*1.96
  }
}
image(longrid[1,],latgrid[1,],temphat,xlab="longtitude",ylab="latitude",col=rainbow(20))
contour(longrid[1,],latgrid[1,],temphat,xlab="longtitude",ylab="latitude")
image(longrid[1,],latgrid[1,],preshat,xlab="longtitude",ylab="latitude",col=rainbow(20))
contour(longrid[1,],latgrid[1,],preshat,xlab="longtitude",ylab="latitude")
image(longrid[1,],latgrid[1,],variance,xlab="longtitude",ylab="latitude",col=rainbow(20))
contour(longrid[1,],latgrid[1,],variance,xlab="longtitude",ylab="latitude")
