# 2017-05-08

data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
save(data,file="MHL4.RData")
bio<-read.table("biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
input<-data[match(rownames(bio),rownames(data)),]

load("MHL2.RData")
pdf("m1.pdf")
boxplot(])
dev.off()
xx<-data[,grep(".6P|X6.P",colnames(data))]
xx<-data[,grep(".7P|X7.P",colnames(data))]
xx<-data[,grep(".|X7.P",colnames(data))]
apply(xx,2,function(x) sum(is.na(x))/nrow(xx))

# automatically select best threshold with 5-fold cross-validation for colon plasma
set.seed(10)
Samping<-CvSampling(75,5)
Lnum<-c()
Mnum<-c()
for(i in 1:5){
  Num<-c()
  for(j in seq(0,1,0.05)){
    counts<-apply(input[,grep(".6P|X6.P",colnames(data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    # counts<-apply(input[,grep(".7P|X7.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>i,na.rm=T)))
    # counts<-apply(input[,grep("NC",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>i,na.rm=T)))
    num<-data.frame(id=j,counts=sum(apply(counts,2,function(x) which.max(x)==2)))
    Num<-rbind(Num,num)
  }
  counts<-apply(input[,grep(".6P|X6.P",colnames(data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>Num[which.max(Num[,2]),1],na.rm=T)))
  counts1<-apply(input[,grep(".6P|X6.P",colnames(data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>Num[which.max(Num[,2]),1],na.rm=T)))
  counts2<-apply(input[,grep(".6P|X6.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>Num[which.max(Num[,2]),1],na.rm=T)))
  Lnum<-c(Lnum,count1=sum(apply(counts1,2,function(x) which.max(x)==2)))
  Mnum<-c(Mnum,count2=sum(apply(counts2,2,function(x) which.max(x)==2)))
}
print(c(Lnum,sum(Lnum)))
print(c(Mnum,sum(Mnum)))

# automatically select best threshold with 5-fold cross-validation for Lung cancer plasma
Samping<-CvSampling(29,5)
Lnum<-c()
Mnum<-c()
for(i in 1:5){
  Num<-c()
  for(j in seq(0,1,0.0001)){
    counts<-apply(input[,grep(".7P|X7.P",colnames(data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    # counts<-apply(input[,grep(".7P|X7.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>i,na.rm=T)))
    # counts<-apply(input[,grep("NC",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>i,na.rm=T)))
    num<-data.frame(id=j,counts=sum(apply(counts,2,function(x) which.max(x)==6)))
    Num<-rbind(Num,num)
  }
  counts1<-apply(input[,grep(".7P|X7.P",colnames(data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>Num[which.max(Num[,2]),1],na.rm=T)))
  counts2<-apply(input[,grep(".7P|X7.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>Num[which.max(Num[,2]),1],na.rm=T)))
  Lnum<-c(Lnum,count1=sum(apply(counts1,2,function(x) which.max(x)==6)))
  Mnum<-c(Mnum,count2=sum(apply(counts2,2,function(x) which.max(x)==6)))
}
print(c(Lnum,sum(Lnum)))
print(c(Mnum,sum(Mnum)/(5*29)))

# automatically select best threshold with 5-fold cross-validation for Normal plasma
Samping<-CvSampling(75,5)
Lnum<-c()
Mnum<-c()
for(i in 1:5){
  Num<-c()
  for(j in seq(0,1,0.005)){
    counts<-apply(input[,grep("NC.P",colnames(data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    num<-data.frame(id=j,counts=sum(apply(counts,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  counts1<-apply(input[,grep("NC.P",colnames(data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>Num[which.max(Num[,2]),1],na.rm=T)))
  counts2<-apply(input[,grep("NC.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>Num[which.max(Num[,2]),1],na.rm=T)))
  Lnum<-c(Lnum,count1=sum(apply(counts1,2,function(x) which.max(x)==10)))
  Mnum<-c(Mnum,count2=sum(apply(counts2,2,function(x) which.max(x)==10)))
}
print(c(Lnum,sum(Lnum)))
print(c(Mnum,sum(Mnum)))

CvSampling<- function(Nobs=30,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}





