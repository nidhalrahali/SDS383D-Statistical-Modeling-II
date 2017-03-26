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

m=100
schoolsum=rep(0,m)
schoolcount=rep(0,m)
for(i in 1:length(mathtest$mathscore)){
  schoolsum[mathtest$school[i]]=schoolsum[mathtest$school[i]]+mathtest$mathscore[i]
  schoolcount[mathtest$school[i]]=schoolcount[mathtest$school[i]]+1
}
t=2000
d=10
eta=10
h=10
theta=gibbssampler(schoolsum,schoolcount,d,eta,h,t)
#use the mean of theta sample to predict the score 
thetamean=colMeans(theta)
plot(thetamean)
predict=mathtest$mathscore
for(i in 1:length(predict)){
  predict[i]=thetamean[mathtest$school[i]]
}
sum((mathtest$mathscore-predict)^2)

schoolmean=schoolsum/schoolcount
kappa=abs(schoolmean-thetamean)/schoolmean
plot(kappa~schoolcount)
