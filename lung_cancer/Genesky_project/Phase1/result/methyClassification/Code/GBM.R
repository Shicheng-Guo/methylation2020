setwd("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/result/methyClassification")
source("class.glasso.r")
source("GSCPackage.R")
library("impute")
library("DAAG")
library("PredictABEL")
library("Deducer")
library("ROCR")
library("MASS")
library("e1071")
library("class")
library("BayesTree")
library("randomForest")
library("gbm")


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

chan<-function(x){
  x[x<0.5]<-0; x[x>=0.5]<-1
  x
}

assess<-function(model,train,test,trainlabel,testlabel){
  t1<-table(chan(predict(model, train)), trainlabel)
  t2<-table(chan(predict(model, test)), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1,na.rm=T)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2,na.rm=T)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}

assess2<-function(model1,model2,trainlabel,testlabel){
  t1<-table(model1,trainlabel)
  t2<-table(model2,testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1,na.rm=T)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2，na.rm=T)
  rlt<-c(sen.train,spe.train,accu.train,sen.test,spe.test,accu.test)
  rlt
}

assess3<-function(model,train,test,trainlabel,testlabel){
  t1<-table((predict(model, train)), trainlabel)
  t2<-table((predict(model, test)), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1,na.rm=T)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2,na.rm=T)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}

assess4<-function(model1,model2,trainlabel,testlabel){
  t1<-table(chan(model1), trainlabel)
  t2<-table(chan(model2), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1,na.rm=T)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2,na.rm=T)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}


data<-read.table("data.txt",head=T,sep="\t",as.is=F,row.names=1)
head(data)
data<-RawNARemove(data)
#dat<-data[,5:19]
dat<-data[,c(2,3,4,5,10,11,12,13,18)]

head(dat)
pheno<-abs(as.numeric(as.factor(data[,20]))-2)
pheno
rlt<-impute.knn(data.matrix(dat),k = 3)
input<-data.frame(rlt$data,pheno)

# glasso cross validation
pheno<-input$pheno
dat<-input[,1:(dim(input)[2]-1)]

# part2
N=10  # times of cross validation
K=5  # fold cross validation
rlt1<-rlt2<-rlt3<-rlt4<-rlt5<-rlt6<-rlt7<-rlt8<-rlt9<-rlt10<-matrix(NA,N*K,6)
j<-1
input<-input[,1:9]


knnrlt<-matrix(NA,9,6)

for(n in 1:N){
  rlt1<-rlt2<-rlt3<-rlt4<-rlt5<-rlt6<-rlt7<-rlt8<-rlt9<-rlt10<-matrix(NA,K,6)
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
    model <- gbm(trainlabel~.,data=train,distribution = "adaboost",n.trees=1500, shrinkage=0.1, bag.fraction=0.5, cv.folds=5)
    pdf("a.pdf")
    best.iter <- gbm.perf(model,method="OOB")
    print(best.iter)
    dev.off()
    model <- gbm(trainlabel~.,data=train,distribution = "adaboost",n.trees=best.iter, shrinkage=0.1, bag.fraction=0.5, cv.folds=5)
    model1 <- predict.gbm(model,train,type="response")
    model2 <- predict.gbm(model,test,type="response")
    rlt5[i,]<-assess2(model1,model2,trainlabel,testlabel)
    
  }
  rltt5[n,]<-colMeans(rlt5,na.rm=T)
}
rltt5
write.table(rltt5,file="random.foreast.accuracy.result.txt",sep="\t",quote=F,row.names=T,col.names=NA)
