gibbssampler=function(ysum,n,t,mu,sigmasq,tausq){
  yvar=sigmasq
  thetavar=yvar*sigmasq*tausq/(yvar+sigmasq*tausq*n)
  updthetamean=mu
  updy=ysum
  for(i in seq(from=1,to=500)){
    updthetamean=(thetavar*updy+yvar*mu)/(yvar+thetavar*n)
    theta=rnorm(1,updthetamean,sqrt(thetavar))
    updy=rnorm(n,theta,sqrt(yvar))
    updy=sum(updy)
    oldthetamean=updthetamean
  }
  theta=matrix(ncol=t)
  for( i in 1:t){
    updthetamean=(thetavar*updy+yvar*mu)/(yvar+thetavar*n)
    thetatemp=rnorm(1,updthetamean,sqrt(thetavar))
    updy=rnorm(n,thetatemp,sqrt(yvar))
    theta[1,i]=thetatemp
    updy=sum(updy)
    oldthetamean=updthetamean
  }
  return=theta
}