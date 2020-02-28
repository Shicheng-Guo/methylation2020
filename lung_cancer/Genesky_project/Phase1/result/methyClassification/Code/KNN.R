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


CvSampling<- function(Nobs,K=5){
  #======================
  # return sample row number: used to sampling the subset to cross-validation predition.
  # rank<-CvSampling(Nobs=dim(data)[1],K=5)
  # Nobs is number of obsevers 
  # train=x[[1]]$train
  # test=x[[1]]$test
  # K=1:10(neibourghs)
  # 
  #====================
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

assess2<-function(model1,model2,trainlabel,testlabel){
  t1<-table(model1,trainlabel)
  t2<-table(model2,testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
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

for(k in 2:10){
  rlt1<-rlt2<-rlt3<-rlt4<-rlt5<-rlt6<-rlt7<-rlt8<-rlt9<-rlt10<-matrix(NA,K,6)
  for(n in 1:N){
    dat<-CvSampling(dim(input)[1],5)
    for(i in 1:K){
      train=input[dat[[i]]$train,]
      trainlabel<-as.factor(pheno[dat[[i]]$train])
      test=input[dat[[i]]$test,]
      testlabel<-pheno[dat[[i]]$test]
      model1<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$train,],cl=pheno[dat[[i]]$train],k=3)
      model2<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$test,],cl=pheno[dat[[i]]$train],k=3)
      rlt2[i,]<-assess2(model1,model2,trainlabel,testlabel)   
    }
  }
  knnrlt[k-1,]<-colMeans(rlt2)
}
write.table(knnrlt,file="knn.accuracy.result.txt",sep="\t",quote=F,row.names=T,col.names=NA)

