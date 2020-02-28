#x: is a data matrix
#y: is a vector, 0 and 1 coded, 1--case, 0--control
#folder: how many folder you want to do the classification

class.glasso <- function(x,y,folder)
{
cv.train <- list()
cv.test <- list()
cv.train.cl <- list()
cv.test.cl <- list()

x.case <- x[y==1,]
x.cntl <- x[y==0,]

y.case <- y[y==1]
y.cntl <- y[y==0]

n.case <- nrow(x.case)
n.cntl <- nrow(x.cntl)

fold1=sample(1:folder,n.case,replace=TRUE, prob=rep(1/folder,folder))
fold2=sample(1:folder,n.cntl,replace=TRUE,prob=rep(1/folder,folder) )

for (i in 1:folder){
  
  train=rbind(x.case[which(fold1!=i),],x.cntl[which(fold2!=i),])  ##sample*features
  test=rbind(x.case[which(fold1==i),],x.cntl[which(fold2==i),])   ##sample*features
  
  train_cl1=matrix(y.case[which(fold1!=i)])
  train_cl2=matrix(y.cntl[which(fold2!=i)])
  
  test_cl1=matrix(y.case[which(fold1==i)])
  test_cl2=matrix(y.cntl[which(fold2==i)])
  
  train_cl=rbind(train_cl1,train_cl2)
  test_cl=rbind(test_cl1,test_cl2)
  
  cv.train[[i]] <- train
  cv.test[[i]] <- test
  cv.train.cl[[i]] <- train_cl
  cv.test.cl[[i]] <- test_cl  
}

library(glmnet)

out <- matrix(0,ncol=6, nrow=folder)

for (k in 1:folder)
{

x <- as.matrix(cv.train[[k]])	
y <- as.numeric(cv.train.cl[[k]])
t_cl <- as.numeric(cv.test.cl[[k]])
test_cl <- t_cl
	
CVD.las <- cv.glmnet(x, y,  family = "binomial", type="class")

tt <- as.matrix(cv.test[[k]])
test <- tt
new <- predict(CVD.las,test, type="class")
new.train <- predict(CVD.las,x, type="class")

o.test <- table(new, test_cl)
o.train <- table(new.train, y)

out[i,1] <- o.test[2,2]/sum(o.test[,2])
out[i,2] <- o.test[1,1]/sum(o.test[,1])
out[i,3] <- (o.test[2,2] + o.test[1,1])/sum(o.test)

out[i,4] <- o.train[2,2]/sum(o.train[,2])
out[i,5] <- o.train[1,1]/sum(o.train[,1])
out[i,6] <- (o.train[2,2] + o.train[1,1])/sum(o.train)
}
colnames(out) <-c("test_sen", "test_spe", "test_acc", "train_sen", "train_spe", "train_acc") 	
n.name <- paste("CV",seq(1,folder), sep="")
rownames(out) <- n.name

return(out)
}