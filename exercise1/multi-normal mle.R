normalmle<-function(x){
  n=dim(x)[1]
  d=dim(x)[2]
  mumle=colSums(x)/n
  sigmamle=matrix(0,nrow=d,ncol=d)
  for(i in 1:n){
    sigmamle=sigmamle+(x[i,]-mumle)%*%t(x[i,]-mumle)
  }
  sigmamle=sigmamle/n
  list(result1=mumle,result2=sigmamle)
}
