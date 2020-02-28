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
library("pROC")


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


data<-read.table("Bigdata.txt",head=T,sep="\t",as.is=F,row.names=1)
head(data)
data<-RawNARemove(data)
dat<-data[,c(2,3,4,5,10,11,12,13,18)]
#dat<-data[,2:19]  # 1) all the methylation site
dim(dat)
head(dat)

# stratified analysis
dat1<-data[which(data$Smoking>0),c(3,4,5,10,15,16,17,18,23,9)]
dat2<-data[which(data$Smoking==0),c(3,4,5,10,15,16,17,18,23,9)]

dat3<-data[which(data$type=="ad"),c(3,4,5,10,15,16,17,18,23,9)]
dat4<-data[which(data$type=="sq"),c(3,4,5,10,15,16,17,18,23,9)]

dat5<-data[which(data$differention=="low"),c(3,4,5,10,15,16,17,18,23,9)]
dat6<-data[c(which(data$differention=="medium"),which(data$differention=="high")),c(3,4,5,10,15,16,17,18,23,9)]

early<-c(which(data$TNM=="IA"),which(data$TNM=="IB"),which(data$TNM=="IIA"),which(data$TNM=="IIB"))
late<-c(which(data$TNM=="IIIA"),which(data$TNM=="IIIB"),which(data$TNM=="IV"))
dat7<-data[early,c(3,4,5,10,15,16,17,18,23,9)]
dat8<-data[late,c(3,4,5,10,15,16,17,18,23,9)]

head(data)

dat<-dat6
{
pheno<-abs(as.numeric(as.factor(dat$pheno))-2)
rlt<-impute.knn(data.matrix(dat),k = 3)
input<-data.frame(rlt$data)
pheno<-input$pheno
dat<-input[,1:(dim(input)[2]-1)]
dim(dat)
head(dat)

# part2
N=100  # times of cross validation
K=5  # fold cross validation
rlt1<-rlt2<-rlt3<-rlt4<-rlt5<-rlt6<-rlt7<-rlt8<-rlt9<-rlt10<-matrix(NA,N*K,6)
r<-1
input<-input[,1:(dim(input)[2]-1)]
head(input)
for(n in 1:N){
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(pheno[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-pheno[dat[[i]]$test]
   
    ### 1 logistic model
    model<-glm(trainlabel~.,family=binomial,data=train)
    predict(model,train,type="response")
    predict(model,test,type="response")
    
    mydata<-train
    mydata$pred<-predict(model, train,type="response")
    mydata$label<-trainlabel
    g <- roc(mydata$label ~ mydata$pred)
    g$auc
    rlt1[r,]<-assess(model,train,test,trainlabel,testlabel,g$acu)
    # 
    r=r+1
  }
}


highpart<-function(rlt1){
  rlt1[match(sort(rowSums(rlt1[,3:6]),decreasing=T)[1:50],rowSums(rlt1[,3:6])),]
}
}

rlt4[,6]<-auc
dim(rlt4)
a<-mean((rlt4)[,6],na.rm=T)
b<-sd((rlt4)[,6],na.rm=T)
n=500

m<-round(a+qnorm(0.95)*b/sqrt(n),2)
n<-round(a,2)
p<-round(a-qnorm(0.95)*b/sqrt(n),2)
paste("Acc=",n,", 95% CI:",p,"-",m,sep="")


low500=rlt4
m500=rlt4
t.test(low500,m500)
m500

rltt<-colMeans(highpart(rlt8),na.rm=T)

rlth<-rbind(rltt1,rltt2,rltt3,rltt6,rltt7,rltt8)
rlth
write.table(rlth,file="accuracy.result.txt",sep="\t",quote=F,row.names=T,col.names=NA)
getwd()


mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
mylogit <- glm(admit ~ gre, data = mydata, family = "binomial")
summary(mylogit)
prob=predict(mylogit,type=c("response"))
mydata$prob=prob
prob

