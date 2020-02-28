# The Second Analysis Pipeline (mix the three different sample together and then find the uniform t and the make the predciton)
# 2017-05-08

setwd("/oasis/tscc/scratch/shg047/monod/hapinfo")
library("ggplot2")
CvSampling<- function(Nobs=30,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}
# data<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.txt",head=T,row.names=1,sep="\t")
#  save(data,file="MHL4.RData")
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
input<-data[match(rownames(bio),rownames(data)),]
# automatically select best threshold with 5-fold cross-validation for colon plasma/lung cancer/normal plasma together.
set.seed(10)
acc1<-c()
acc2<-c()
acc3<-c()
for(loop in 1:10){
Samping<-CvSampling(29,5)
Lnum1<-c()
Lnum2<-c()
Lnum3<-c()
for(i in 1:5){
  Num<-c()
  for(j in seq(0,1,0.005)){
    counts1<-apply(input[,grep(".6P|X6.P",colnames(data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    counts2<-apply(input[,grep(".7P|X7.P",colnames(data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    counts3<-apply(input[,grep("NC.P",colnames(data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    num<-data.frame(id=j,
                    c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                    c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                    c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  best<-Num[which.max(rowSums(Num[,2:4])),1]
  countm1<-apply(input[,grep(".6P|X6.P",colnames(data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm2<-apply(input[,grep(".7P|X7.P",colnames(data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm3<-apply(input[,grep("NC.P",colnames(data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  Lnum1<-c(Lnum1,c=sum(apply(countm1,2,function(x) which.max(x)==2)))
  Lnum2<-c(Lnum2,c=sum(apply(countm2,2,function(x) which.max(x)==6)))
  Lnum3<-c(Lnum3,c=sum(apply(countm3,2,function(x) which.max(x)==10)))
  
  acc1<-rbind(acc1,c(best,sum(apply(countm1,2,function(x) which.max(x)==2))/(length(Samping[[i]]$test))))
  acc2<-rbind(acc2,c(best,sum(apply(countm2,2,function(x) which.max(x)==6))/(length(Samping[[i]]$test))))
  acc3<-rbind(acc3,c(best,sum(apply(countm3,2,function(x) which.max(x)==10))/(length(Samping[[i]]$test))))
}
print(loop)
}
save.image()
senthresplot<-function(){
acc<-data.frame(x=acc[,1],y=acc[,2])
acc$bin<- cut(acc[,1], c(seq(0,1,0.01)))
ggplot(acc) + geom_boxplot(aes(bin, y))
}
senthresplot(acc1)
ggsave("colon-threshold-acc.pdf")
senthresplot(acc2)
ggsave("lung-threshold-acc.pdf")
senthresplot(acc3)
ggsave("normal-threshold-acc.pdf")
