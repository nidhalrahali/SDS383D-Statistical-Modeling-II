mathtest <- read.csv("~/GitHub/SDS383D-course-work/exercise4/mathtest.csv")
boxplot(mathscore~school,data=mathtest)

# fit a model for mathscore versus school
lm1 = lm(mathscore ~ school, data=mathtest)
pred = predict(lm1, mathtest)
plot(pred, mathtest$mathscore)
sum( (mathtest$mathscore - pred)^2 )

#compare with the grandmean of all scores
grandmean = mean(mathtest$mathscore)
sum( (mathtest$mathscore - grandmean)^2 )
