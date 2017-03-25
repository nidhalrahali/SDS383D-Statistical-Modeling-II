mathtest <- read.csv("~/GitHub/SDS383D-course-work/exercise4/mathtest.csv")
boxplot(mathscore~school,data=mathtest)
summary(mathtest)
# fit a simple linear model for mathscore versus school
lm1 = lm(mathscore ~ school, data=mathtest)
pred = predict(lm1, mathtest)
plot(pred, mathtest$mathscore)
sum( (mathtest$mathscore - pred)^2 )

#compare with the grandmean of all scores
grandmean = mean(mathtest$mathscore)
sum( (mathtest$mathscore - grandmean)^2 )

schoolscore=rep(0,max(mathtest$school))
schoolsize=rep(0,max(mathtest$school))
for(i in 1:length(mathtest$school)){
  schoolscore[mathtest$school[i]]=schoolscore[mathtest$school[i]]+mathtest$mathscore[i]
  schoolsize[mathtest$school[i]]=schoolsize[mathtest$school[i]]+1
}
t=1000
mu=grandmean
sigmasq=1
tausq=1
theta=matrix(nrow=max(mathtest$school),ncol=t)
for(i in 1:length(schoolscore)){
theta[i,]=gibbssampler(totalschoolscore[i],schoolsize[i],t,mu,sigmasq,tausq)
}
thetamean=rowMeans(theta)
plot(thetamean)
predict=mathtest$mathscore
for(i in 1:length(predict)){
  predict[i]=thetamean[mathtest$school[i]]
}
sum((mathtest$mathscore-predict)^2)
