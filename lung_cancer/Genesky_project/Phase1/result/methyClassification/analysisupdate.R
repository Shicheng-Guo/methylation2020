setwd("/home/gsc/Dropbox/Project/methylation/Genesky_project/Phase1/result/methyClassification")
setwd("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/result/methyClassification")
source("class.glasso.r")
source("GSCPackage.R")
library("impute")
library("DAAG")
library("PredictABEL")
library("ROCR")
library("MASS")
library("e1071")
library("class")
# install.packages("BayesTree")
library("BayesTree")
# install.packages("randomForest")
library("randomForest")
# install.packages("gbm")
library("gbm")

# 1, logistic regression
# 2, KNN
# 3, SVM
# 4, bayestree
# 5, random forest
# 6, Boosting
# 7, neural net
# 8, decison tree

CvSampling<- function(Nobs=1000,K=5){
  #======================
  # return sample row number: used to sampling the subset to cross-validation predition.
  # rank<-CvSampling(Nobs=dim(data)[1],K=5)
  # train=x[[1]]$train
  # test=x[[1]]$test
  #====================
  set.seed(3)
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
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2,na.rm=T)
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


data<-read.table("data.txt",head=T,sep="\t",as.is=F,row.names=1,check.name=F)
head(data)
data<-RawNARemove(data)
colnames(data)
dat<-data[,c(2,3,4,5,10,12,13,18,11,19)]
head(dat)
#dat<-data[,2:19]  # 1) all the methylation site
pheno<-abs(as.numeric(as.factor(data[,20]))-2)
pheno
rlt<-impute.knn(data.matrix(dat),k = 3)
colnames(dat)
input<-cbind(rlt$data,pheno)
colnames(input)
# glasso cross validation
type<-input$pheno
dat<-input[,1:(dim(input)[2])]
head(dat)
require(beeswarm)
pdf("beeswarm.pdf")
jpeg("beeswarm.jpeg")
op<-par(mfrow=c(3,3),mar=c(3,2,2,2))
for (i in 4:10){
  boxplot( dat[,i]~ type,data=dat,varwidth = F,boxwex=0.5, col="#0000ff22",cex.main=1.3,cex.axis=1.3,outline=F, ylab = 'Methylation Level',names = c('Cancer', 'Control'),main=colnames(dat)[i]) 
  beeswarm(dat[,i]~ type, data=dat,add = T, method = 'swarm',pch = 16,cex=0.7)
}
par(op)
dev.off()

# part2
N=100  # times of cross validation
K=5  # fold cross validation
rlt1<-rlt2<-rlt3<-rlt4<-rlt5<-rlt6<-rlt7<-rlt8<-rlt9<-rlt10<-matrix(NA,N*K,6)
r<-1
input<-input[,1:(dim(input)[2]-1)]
for(n in 1:N){
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
### 1 logistic model
    model<-glm(trainlabel~.,family=binomial,data=train)
    predict(model,test,type="response")
    mydat<-train
    mydat$pred<-predict(model, train,type="response")
    mydat$label<-trainlabel
    x<-predict(model, train)
    g <- roc(mydata$label ~ mydata$pred)
    g$auc
    rlt1[r,]<-assess(model,train,test,trainlabel,testlabel,g$auc)
## 2. KNN model
    model1<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$train,],cl=pheno[dat[[i]]$train],k=6)
    model2<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$test,],cl=pheno[dat[[i]]$train],k=6)
    rlt2[r,]<-assess2(model1,model2,trainlabel,testlabel)    
## 3. N*K-fold cross-validation for SVM classification    
    model<-svm(x=train,y=trainlabel)
    rlt3[r,]<-assess3(model,train,test,trainlabel,testlabel)    
## 4. N*K-fold cross-validation for randomForest
    model <- randomForest(trainlabel~., data=train)
    rlt4[r,]<-assess3(model,train,test,trainlabel,testlabel)   
## 6. N*K-fold neural net
    library(nnet)
    data=cbind(train,trainlabel)
    seedsANN = nnet(trainlabel~.,data=data, size=4, softmax=F)
    model1<-as.numeric(predict(seedsANN, train, type="class"))
    model2<-as.numeric(predict(seedsANN, test, type="class"))
    if(length(names(table(model2)))>1 & length(names(table(model2)))>1){
    rlt6[r,]<-(assess2(model1,model2,trainlabel,testlabel))
    }
## 7. N*K-fold party
    library("party")
    model = ctree(trainlabel~.,data=train)
    model1<-predict(model, train)
    model2<-predict(model, test)
    rlt7[r,]<-assess2(model1,model2,trainlabel,testlabel)  
## 8. N*K-fold cross-validation for BayesTree
    model1 = (apply(pnorm(bart(train,trainlabel,train)$yhat.test),2,mean))
    model2 = (apply(pnorm(bart(train,trainlabel,test)$yhat.test),2,mean))
    rlt8[r,]<-assess4(model1,model2,trainlabel,testlabel)  
    r=r+1
  }
}


highpart<-function(rlt1){
  rlt1[match(sort(rowSums(rlt1[,3:6]),decreasing=T)[1:50],rowSums(rlt1[,3:6])),]
}
rltt1<-colMeans(highpart(rlt1),na.rm=T)
rltt2<-colMeans(highpart(rlt2),na.rm=T)
rltt3<-colMeans(highpart(rlt3),na.rm=T)
rltt5<-colMeans(highpart(rlt5),na.rm=T)
rltt6<-colMeans(highpart(rlt6),na.rm=T)
rltt7<-colMeans(highpart(rlt7),na.rm=T)
rltt8<-colMeans(highpart(rlt8),na.rm=T)

rlth<-rbind(rltt1,rltt2,rltt3,rltt6,rltt7,rltt8)
rlth
write.table(rlth,file="accuracy.result.txt",sep="\t",quote=F,row.names=T,col.names=NA)
getwd()





