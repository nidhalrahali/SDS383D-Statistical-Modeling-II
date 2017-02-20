source('~/GitHub/SDS383D-course-work/exercise3/weight.R')

testmodel=function(ytrain,xtrain,ytest,xtest,h,t){
  #normalize all data
  ymean=mean(ytrain)
  ynorm=ytrain-matrix(ymean,ncol=length(ytrain))
  xnorm=xtrain-matrix(mean(xtrain),ncol=length(xtrain))
  xtestnorm=xtest-matrix(mean(xtrain),ncol=length(xtest))
  
  # fit smoother and record the error
  error=matrix(ncol=length(ytrain))
  for(i in 1:length(ytrain)){
    error[i]=ynorm %*% t(weight(xtestnorm[i],xnorm,h,type=t))+ymean-ytest[i]
  }
  return(error)
}