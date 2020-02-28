




saminfo<-read.table("/home/shg047/oasis/monod/saminfo.txt",sep="\t")
data<-read.table("/home/shg047/oasis/monod/hapinfo/monod.mhl.may5.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(data)

colnames(data)[grep("STL",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("STL",colnames(data))],saminfo[,1]),2])
colnames(data)

x1<-grep("Colon|CTT|6-P|6-T",colnames(data), ignore.case=T)
x2<-grep("Lung|7-P|7-T",colnames(data), ignore.case=T)
x3<-grep("Pancr|PC-P|PC-T",colnames(data), ignore.case=T)
x4<-grep("NC-|WB",colnames(data), ignore.case=T)

colnames(data)[x1]
colnames(data)[x2]
colnames(data)[x3]
colnames(data)[x4]

data<-data[,c(x1,x2,x3,x4)]

colnames(data)[grep("6-P-",colnames(data))]<-"6-P";
colnames(data)[grep("7-P-",colnames(data))]<-"7-P";
colnames(data)[grep("6-T-",colnames(data))]<-"6-T";
colnames(data)[grep("7-T-",colnames(data))]<-"7-T";
colnames(data)[grep("PC-P-",colnames(data))]<-"PC-P";
colnames(data)[grep("PC-T-",colnames(data))]<-"PC-T";
colnames(data)[grep("NC-P-",colnames(data))]<-"NC-P";
colnames(data)
save(data,file="Tissue-specific-MHBs-Figure5-may2016.RData")

library("impute")
f2<-RawNARemove(data,missratio=0.3)
data<-impute.knn(data.matrix(f2))$data


source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("stringr")
biocLite("randomForest")
biocLite("impute")
biocLite("rpart")
biocLite("e1071")
biocLite("biclust")

library("Biostrings")
library("stringr")
library("randomForest")
library("impute")
library("rpart")
library("e1071")
library("biclust")

setwd("/home/shg047/monod/mhl")
sam<-read.table("../sampleinfo.sort.txt",sep="\t",as.is=T)
h1<-read.table("N37_WB_WGBS_RRBS_dRRBS_merged__mld_blocks_stringent_mhl_matrix.txt",sep="\t",head=T,as.is=T,check.names=F,row.names=1)
h2<-read.table("N37_WB_WGBS_SeqCap_BSPP_merged__mld_blocks_stringent_mhl_matrix.txt",sep="\t",head=T,as.is=T,check.names=F,row.names=1)
phen1<-sam[match(colnames(h1),sam[,1]),2]
phen2<-sam[match(colnames(h2),sam[,1]),2]
# check the significant differential regions

PairWilPValue<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],5)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(wilcox.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=F, na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-mean(as.numeric(data[i,x1]),na.rm=T)-mean(as.numeric(data[i,x2]),na.rm=T)
      output[i,3]<-"wilcox"
      output[i,4]<-mean(as.numeric(data[i,x1]),na.rm=T)
      output[i,5]<-mean(as.numeric(data[i,x2]),na.rm=T)
      print(i)
    }
  }
  out<-cbind(rownames(data),output)
  out
}

x1<-which(phen1==names(table(phen1))[1])   # type 1, cancer or sensitive
x2<-which(phen1==names(table(phen1))[2])   # type 2, normal or resistant
rlt1<-PairWilPValue(h1,x1,x2)
x1<-which(phen2==names(table(phen2))[1])   # type 1, cancer or sensitive
x2<-which(phen2==names(table(phen2))[2])   # type 2, normal or resistant
rlt2<-PairWilPValue(h2,x1,x2)
rlt1<-data.frame(rlt1)
rlt2<-data.frame(rlt2)
colnames(rlt1)<-c("Bin","P","Delta","test","cancer","normal")
colnames(rlt2)<-c("Bin","P","Delta","test","cancer","normal")
rlt1$P<-as.numeric(as.character(rlt1$P))
rlt2$P<-as.numeric(as.character(rlt2$P))
rlt1$Delta<-as.numeric(as.character(rlt1$Delta))
rlt2$Delta<-as.numeric(as.character(rlt2$Delta))
rlt1$normal<-as.numeric(as.character(rlt1$normal))
rlt2$normal<-as.numeric(as.character(rlt2$normal))

rlt3<-subset(rlt1,P<0.05/nrow(rlt1) & Delta>0 & normal<0.1) # 34
rlt4<-subset(rlt2,P<0.05/nrow(rlt2)) # 2

rlt3<-subset(rlt1,P<0.05/nrow(rlt1) & Delta>0 & normal<0.2) # 36
rlt4<-subset(rlt2,P<0.05/nrow(rlt2) & Delta>0 & normal<0.2) # 8

rlt5<-subset(rlt1,P<0.05/nrow(rlt1) & Delta>0 & normal<0.3) # 36
rlt6<-subset(rlt2,P<0.05/nrow(rlt2) & Delta>0 & normal<0.3) # 13
dim(rlt5)
dim(rlt6)

write.table(rlt1,file="N37_WB_WGBS_RRBS_dRRBS_merged__mld_blocks_stringent_mhl_matrix.txt.Pvalue.txt",row.names=T,col.names=NA,sep="\t",quote=F)
write.table(rlt2,file="N37_WB_WGBS_SeqCap_BSPP_merged__mld_blocks_stringent_mhl_matrix.txt.Pvalue.txt",row.names=T,col.names=NA,sep="\t",quote=F)
write.table(rlt3,file="N37_WB_WGBS_RRBS_dRRBS_merged__mld_blocks_stringent_mhl_matrix.txt.sig.txt",row.names=T,col.names=NA,sep="\t",quote=F)
write.table(rlt4,file="N37_WB_WGBS_SeqCap_BSPP_merged__mld_blocks_stringent_mhl_matrix.txt.sig.txt",row.names=T,col.names=NA,sep="\t",quote=F)


# check the significant differential regions
Cvsampling<- function(phen,k=2){
  #======================
  # suppose phen has two type: used to sampling the subset to cross-validation predition with same proportion of case and control in train and test
  # rank<-CvSampling(nobs=dim(data)[1],k=5)
  # train=x[[1]]$train
  # test=x[[1]]$test
  #====================
  rid<-list()
  ll1<-list()
  ll2<-list()
  l1<-which(phen==names(table(phen)[1]))
  l1<-l1[order(runif(length(l1)))]
  l2<-which(phen==names(table(phen)[2]))
  l2<-l2[order(runif(length(l2)))]
  kk1 <- as.integer(length(l1)*seq(1,k-1)/k)  # sample size in each node
  kk2 <- as.integer(length(l2)*seq(1,k-1)/k)  # sample size in each node
  kk1 <- matrix(c(0,rep(kk1,each=2),length(l1)),ncol=2,byrow=TRUE)
  kk2 <- matrix(c(0,rep(kk2,each=2),table(phen)[2]),ncol=2,byrow=TRUE)
  kk1[,1] <- kk1[,1]+1
  kk2[,1] <- kk2[,1]+1
  ll1 <- lapply(seq.int(k),function(x,kk=kk1,d=l1) list(train=d[!(seq(d) %in% seq(kk[x,1],kk[x,2]))], test=d[seq(kk[x,1],kk[x,2])]))
  ll2 <- lapply(seq.int(k),function(x,kk=kk2,d=l2) list(train=d[!(seq(d) %in% seq(kk[x,1],kk[x,2]))], test=d[seq(kk[x,1],kk[x,2])]))
  rid<-lapply(seq.int(k),function(x,y=ll1,z=ll2) list(train=c(y[[x]]$train,z[[x]]$train),test=c(y[[x]]$test,z[[x]]$test)))
  return(rid)
}

setwd("/home/shg047/monod/mhl")
sam<-read.table("../sampleinfo.sort.txt",sep="\t",as.is=T)
h1<-read.table("N37_WB_WGBS_RRBS_dRRBS_merged__mld_blocks_stringent_mhl_matrix.txt",sep="\t",head=T,as.is=T,check.names=F,row.names=1)
h2<-read.table("N37_WB_WGBS_SeqCap_BSPP_merged__mld_blocks_stringent_mhl_matrix.txt",sep="\t",head=T,as.is=T,check.names=F,row.names=1)

# check the raw data of mhl for the significant differential mhl
sig1<-read.table(file="N37_WB_WGBS_RRBS_dRRBS_merged__mld_blocks_stringent_mhl_matrix.txt.sig.txt",head=T,row.names=1,sep="\t",as.is=T)
sig2<-read.table(file="N37_WB_WGBS_SeqCap_BSPP_merged__mld_blocks_stringent_mhl_matrix.txt.sig.txt",head=T,row.names=1,sep="\t",as.is=T)

u1<-h1[match(sig1[,1],rownames(h1)),]
u2<-h2[match(sig2[,1],rownames(h2)),]
write.table(u1,file="N37_WB_WGBS_RRBS_dRRBS_merged__mld_blocks_stringent_mhl_matrix.txt.sig.raw.mhl.txt",row.names=T,col.names=NA,sep="\t",quote=F)
write.table(u2,file="N37_WB_WGBS_SeqCap_BSPP_merged__mld_blocks_stringent_mhl_matrix.txt.sig.raw.mhl.txt",row.names=T,col.names=NA,sep="\t",quote=F)

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}  

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
} 

nu1<-ColNARemove(u1)
nu1<-RawNARemove(nu1)
nu1<-data.matrix(nu1)
nnu1<-impute.knn(nu1)$data


data<-t(nnu1)
phen1<-sam[match(colnames(nnu1),sam[,1]),2]
pheno=phen1

nu1<-ColNARemove(u2)
nu1<-RawNARemove(nu1)
nu1<-data.matrix(nu1)
nnu1<-impute.knn(nu1)$data
data<-t(nnu1)
phen1<-sam[match(colnames(nnu1),sam[,1]),2]
pheno=phen1

pheno=y
data=data
N=10  # times of cross validation
K=2  # fold cross validation
rlt1<-matrix(NA,N*K,6)
rlt2<-matrix(NA,N*K,6)
rlt3<-matrix(NA,N*K,6)
j<-1

for(n in 1:N){
  index<-Cvsampling(pheno,k=K)
  for(i in 1:K){
    train=data.frame(data[index[[i]]$train,])
    trainlab=pheno[index[[i]]$train]
    test=data.frame(data[index[[i]]$test,])
    testlab=pheno[index[[i]]$test]
    
    # for randomeForests
    model <- randomForest(as.factor(trainlab)~., data=train,importance=T,na.action="na.omit") # as.factor(y)
    predict1 <- predict(model, train)
    predict2 <- predict(model, test)
    t1<-table(predict1, trainlab)
    t2<-table(predict2, testlab)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt1[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    
    # for rpart
    model <- rpart(as.factor(trainlab) ~.,data=train,method="class")
    predict1 <- predict(model, train,type="class")
    predict2 <- predict(model, test, type="class")
    t1<-table(predict1, trainlab)
    t2<-table(predict2, testlab)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt2[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    
    # for svm
    model <- svm(as.factor(trainlab) ~.,data=train)
    predict1 <- predict(model, train)
    predict2 <- predict(model, test)
    t1<-table(predict1, trainlab)
    t2<-table(predict2, testlab)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt3[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    
    # for LR
    print(j)
    j=j+1
  }
}

rrlt1<-colSums(rlt1)/(K*N)
rrlt2<-colSums(rlt2)/(K*N)
rrlt3<-colSums(rlt3)/(K*N)

out1<-rbind(rrlt1,rrlt2,rrlt3)
out1<-round(out1,3)
colnames(out1)<-paste(rep(c("specificity","sensitivity","accuray"),2),rep(c("train","test"),each=3),sep="_")
rownames(out1)<-c("random forest","recursive partitioning trees","support vector machine")
write.table(out1,file="sen.spe.acc.file1.txt",sep="\t",col.names=NA,row.names=T,quote=F)


out2<-rbind(rrlt1,rrlt2,rrlt3)
out2<-round(out2,3)
colnames(out2)<-paste(rep(c("specificity","sensitivity","accuray"),2),rep(c("train","test"),each=3),sep="_")
rownames(out2)<-c("random forest","recursive partitioning trees","support vector machine")
write.table(out2,file="sen.spe.acc.file1.txt",sep="\t",col.names=NA,row.names=T,quote=F)




library("randomForest")
pred2<-mydata
lab2<-as.factor(colnames(mydata))
bestmtry<-tuneRF(x=t(pred2),y=lab2,ntreeTry=1000, stepFactor=2, improve=0.05,trace=TRUE, plot=TRUE, doBest=T)
best.mtry=bestmtry$mtry
rf.model1<-randomForest(x=t(pred2),y=lab2,mtry=best.mtry)
rf.model1
rf.importance<-rf.importance+rf.model1$importance/100

















