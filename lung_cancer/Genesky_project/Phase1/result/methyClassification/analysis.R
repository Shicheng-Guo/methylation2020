setwd("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/result/methyClassification")
# 1, logistic regression
# 2, KNN
# 3, SVM
# 4, bayestree
# 5, random forest
# 6, Boosting
# 7, Fuzzy Rule-based Systems
# source("class.glasso.r")
source("GSCPackage.R")
library("impute")
library("DAAG")
library("PredictABEL")
library("Deducer")
library("Rcmdr")
library("ROCR")
library("MASS")

data<-read.table("data.txt",head=T,sep="\t",as.is=F,row.names=1)
head(data)
data<-RawNARemove(data)
dat<-data[,5:19]

head(dat)
pheno<-abs(as.numeric(as.factor(data[,20]))-2)
pheno
rlt<-impute.knn(data.matrix(dat),k = 3)
input<-data.frame(rlt$data,pheno)

library("class")
CvSampling<- function(Nobs=1000,K=5){
  #======================
  # return sample row number: used to sampling the subset to cross-validation predition.
  # rank<-CvSampling(Nobs=dim(data)[1],K=5)
  # train=x[[1]]$train
  # test=x[[1]]$test
  #====================
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

N=100  # times of cross validation
K=5  # fold cross validation
rlt<-matrix(NA,N*K,6)
j<-1
input<-input[,1:14]
for(n in 1:N){
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    dat<-CvSampling(dim(input)[1],5)
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
    
    # 1. logistic regressopm   
    model<-glm(trainlabel~.,family=binomial,data=train)
    predict(model,train,type="response")
    predict(model,test,type="response")
    predict(model, train)
    x<-predict(model, train)
    chan<-function(x){
      x[x<0.5]<-0; x[x>=0.5]<-1
      x
    }
    t1<-table(chan(predict(model, train)), trainlabel)
    t2<-table(chan(predict(model, test)), testlabel)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    j=j+1
    print(j)
  }
}
colMeans(rlt)

# 2. KNN
N=1000  # times of cross validation
K=5  # fold cross validation
rlt<-matrix(NA,N*K,6)
j<-1
input<-input[,1:14]
for(n in 1:N){
dat<-CvSampling(dim(input)[1],5)
for(i in 1:K){
  xknn1<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$train,],cl=pheno[dat[[i]]$train],k=6)
  xknn2<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$test,],cl=pheno[dat[[i]]$train],k=6)
  t1<-table(xknn1,pheno[dat[[i]]$train])
  t2<-table(xknn2,pheno[dat[[i]]$test])
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  j=j+1
  print(j)
} 
}
colMeans(rlt)


#########
## 3. N*K-fold cross-validation for SVM classification
#########
library(e1071)
N=1000  # times of cross validation
K=5  # fold cross validation
rlt<-matrix(NA,N*K,6)
j<-1
input<-input[,1:14]
for(n in 1:N){
  order<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
    model<-svm(x=train,y=trainlabel)
    t1<-table(predict(model, train), trainlabel)
    t2<-table(predict(model, test), testlabel)
    t1
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    j=j+1
    print(j)
  }
}
colMeans(rlt)




#########
## 4. N*K-fold cross-validation for BayesTree
#########
library("BayesTree")
library(verification)
N=10  # times of cross validation
K=5  # fold cross validation
rlt<-matrix(NA,N*K,6)
j<-1
input<-input[,1:14]
for(n in 1:N){
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
    chan<-function(x){
      x[x<0.5]<-0; x[x>=0.5]<-1
      x
    }
    predict1 = chan(apply(pnorm(bart(train,trainlabel,train)$yhat.test),2,mean))
    predict2 = chan(apply(pnorm(bart(train,trainlabel,test)$yhat.test),2,mean))
    t1<-table(predict1, trainlabel)
    t2<-table(predict2, testlabel)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    j=j+1
    print(j)
  }
}
colMeans(rlt)



#########
## 5. N*K-fold cross-validation for randomForest
#########
# 
library("randomForest")
N=10  # times of cross validation
K=5  # fold cross validation
rlt<-matrix(NA,N*K,6)
j<-1
input<-input[,1:14]
for(n in 1:N){
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
    chan<-function(x){
      x[x<0.5]<-0; x[x>=0.5]<-1
      x
    }
    model <- randomForest(trainlabel~., data=train)
    predict1 <- predict(model, train)
    predict2 <- predict(model, test)
    t1<-table(predict1, trainlabel)
    t2<-table(predict2, testlabel)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    j=j+1
    print(j)
  }
}
colMeans(rlt)



#6.Boosting
library("gbm")
N=10  # times of cross validation
K=5  # fold cross validation
rlt<-matrix(NA,N*K,6)
j<-1
input<-input[,1:14]
for(n in 1:N){
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
    chan<-function(x){
      x[x<0.5]<-0; x[x>=0.5]<-1
      x
    }
    model <- gbm(trainlabel~.,data=train,distribution = "bernoulli")
    predict1 <- predict.gbm(model,train,100,type="response")
    predict1
    predict2 <- predict.gbm(model,test,100)
    t1<-table(predict1, trainlabel)
    t2<-table(predict2, testlabel)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    j=j+1
    print(j)
  }
}
colMeans(rlt)

##7.Fuzzy Rule-based Systems


# load(input)  input was creat from server since in my local computer impute can't installed
setwd("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/result/methyClassification")
rm(list=ls())
load("input.RDdata")
ls()

# glasso cross validation
pheno<-input$pheno
dat<-input[,1:(dim(input)[2]-1)]
out<-class.glasso(dat,pheno,5)
out




