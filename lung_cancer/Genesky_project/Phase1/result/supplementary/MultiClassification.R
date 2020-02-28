
#######  multiclass classification based on Generalized Boosted Models  ##########
library(gbm)
library(caret)
data(iris)
fitControl <- trainControl(method="repeatedcv",number=5,repeats=1,verboseIter=TRUE)
set.seed(825)
gbmFit <- train(Species ~ ., data=iris,method="gbm",trControl=fitControl,verbose=FALSE)
gbmFit


#######  multiclass classification based on Support Vector Machine  ##########
########## http://rss.acs.unt.edu/Rdoc/library/kernlab/html/ksvm.html #####3
data(iris)
## Create a kernel function using the build in rbfdot function
rbf <- rbfdot(sigma=0.1)
rbf
## train a bound constraint support vector machine
irismodel <- ksvm(Species~.,data=iris,type="C-bsvc",kernel=rbf,C=10,prob.model=TRUE)
irismodel
## get fitted values
fitted(irismodel)
## Test on the training set with probabilities as output
predict(irismodel, iris[,-5], type="probabilities")
## Demo of the plot function
x <- rbind(matrix(rnorm(120),,2),matrix(rnorm(120,mean=3),,2))
y <- matrix(c(rep(1,60),rep(-1,60)))
svp <- ksvm(x,y,type="C-svc")
plot(svp,data=x)


#######  mvrpart  ##########







