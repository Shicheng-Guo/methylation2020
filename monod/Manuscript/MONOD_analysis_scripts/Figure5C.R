# 2017-05-08
# The Second Analysis Pipeline (leave 5 samples for each group and then mix the three different sample together and then find the uniform t and the make the predciton)

setwd("/oasis/tscc/scratch/shg047/monod/hapinfo")
library("ggplot2")

CvSampling<- function(Nobs=29,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

################################################################
################  Model Building ##############################
################################################################
data<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.txt",head=T,row.names=1,sep="\t")
save(data,file="MHL4.RData")
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
  
# Response: Random select 5 samples for each group, leave them alone and then apply our modle to these samples.
out1<-grep(".6P|X6.P",colnames(data))[sample(1:30,5)]
out2<-grep(".7P|X7.P",colnames(data))[sample(1:30,5)]
out3<-grep("NC.P",colnames(data))[sample(1:75,5)]
test<-data[,c(out1,out2,out3)]
test<-test[match(rownames(bio),rownames(test)),]

# Cross-validation process
Data<-data[,-c(out1,out2,out3)]
# automatically select best threshold with 5-fold cross-validation for colon plasma/lung cancer/normal plasma together.
input<-Data[match(rownames(bio),rownames(Data)),]
set.seed(51)
k=2                    # split the sample to two parts, one is for train and the other is for test. 
acc1<-c()
acc2<-c()
acc3<-c()
Best<-c()
Samping<-CvSampling(24,k)
Lnum1<-c()
Lnum2<-c()
Lnum3<-c()
for(i in 1:k){
  Num<-c()
  for(j in seq(0,1,0.001)){
    counts1<-apply(input[,grep(".6P|X6.P",colnames(Data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    counts2<-apply(input[,grep(".7P|X7.P",colnames(Data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    counts3<-apply(input[,grep("NC.P",colnames(Data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    num<-data.frame(id=j,
                      c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                      c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                      c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  best<-Num[which.max(rowSums(Num[,2:4])),1]
  countm1<-apply(input[,grep(".6P|X6.P",colnames(Data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm2<-apply(input[,grep(".7P|X7.P",colnames(Data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm3<-apply(input[,grep("NC.P",colnames(Data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  Lnum1<-c(Lnum1,c=sum(apply(countm1,2,function(x) which.max(x)==2)))
  Lnum2<-c(Lnum2,c=sum(apply(countm2,2,function(x) which.max(x)==6)))
  Lnum3<-c(Lnum3,c=sum(apply(countm3,2,function(x) which.max(x)==10)))
    
  acc1<-rbind(acc1,c(best,sum(apply(countm1,2,function(x) which.max(x)==2))/(length(Samping[[i]]$test))))
  acc2<-rbind(acc2,c(best,sum(apply(countm2,2,function(x) which.max(x)==6))/(length(Samping[[i]]$test))))
  acc3<-rbind(acc3,c(best,sum(apply(countm3,2,function(x) which.max(x)==10))/(length(Samping[[i]]$test))))
  
  Best<-c(Best,best)
}
Best
cc1<-apply(test[,grep(".6P|X6.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
cc2<-apply(test[,grep(".7P|X7.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
cc3<-apply(test[,grep("NC.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
ccc1=sum(apply(cc1,2,function(x) which.max(x)==2))/5
ccc2=sum(apply(cc2,2,function(x) which.max(x)==6))/5
ccc3=sum(apply(cc3,2,function(x) which.max(x)==10))/5
ccc1
ccc2
ccc3
######################################################################
################  Model Stability ####################################
######################################################################
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
starlt<-stability(data,bio,rep=1000)
save(starlt,file="starlt.RData")
# starlt
#        [,1] [,2] [,3] [,4]
# rlt 0.3605  0.8  0.8  0.8
# rlt 0.3390  0.8  1.0  1.0
load("starlt.RData")
senthresplot<-function(starlt){
  acc<-data.frame(x=acc[,1],ACC=acc[,2])
  acc$MHL<- cut(acc[,1], c(seq(0,1,0.1)))
  p <-ggplot(acc,aes(factor(MHL),ACC)) + 
    geom_boxplot(aes(MHL, ACC))+
    ylim(0,1)+
    geom_point(position = position_jitter(width = 0.2))+
                 theme_bw()+
                 theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank())    
}

acc=starlt[,c(1,2)]
senthresplot(acc)
ggsave("colon-threshold-acc-2.pdf")
acc=starlt[,c(1,3)]
senthresplot(acc)
ggsave("lung-threshold-acc-2.pdf")

acc=starlt[,c(1,4)]
senthresplot(acc)
ggsave("normal-threshold-acc-2.pdf")

stability<-function(data,bio,rep){
  Rlt<-c()
  for(loop in 1:rep){
  print(loop)
  out1<-grep(".6P|X6.P",colnames(data))[sample(1:30,5)]
  out2<-grep(".7P|X7.P",colnames(data))[sample(1:29,5)]
  out3<-grep("NC.P",colnames(data))[sample(1:75,5)]
  test<-data[,c(out1,out2,out3)]
  test<-test[match(rownames(bio),rownames(test)),]
  Data<-data[,-c(out1,out2,out3)]
  # automatically select best threshold with 5-fold cross-validation for colon plasma/lung cancer/normal plasma together.
  input<-Data[match(rownames(bio),rownames(Data)),]
  k=2
  acc1<-c()
  acc2<-c()
  acc3<-c()
  Best<-c()
  Samping<-CvSampling(24,k)
  Lnum1<-c()
  Lnum2<-c()
  Lnum3<-c()
  for(i in 1:k){
    Num<-c()
    for(j in seq(0,1,0.05)){
      counts1<-apply(input[,grep(".6P|X6.P",colnames(input))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
      counts2<-apply(input[,grep(".7P|X7.P",colnames(input))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
      counts3<-apply(input[,grep("NC.P",colnames(input))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
      num<-data.frame(id=j,
                      c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                      c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                      c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
  Num<-rbind(Num,num)
  }
  best<-Num[which.max(rowSums(Num[,2:4])),1]
  countm1<-apply(input[,grep(".6P|X6.P",colnames(input))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm2<-apply(input[,grep(".7P|X7.P",colnames(input))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm3<-apply(input[,grep("NC.P",colnames(input))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  Lnum1<-c(Lnum1,c=sum(apply(countm1,2,function(x) which.max(x)==2)))
  Lnum2<-c(Lnum2,c=sum(apply(countm2,2,function(x) which.max(x)==6)))
  Lnum3<-c(Lnum3,c=sum(apply(countm3,2,function(x) which.max(x)==10)))
    
  acc1<-rbind(acc1,c(best,sum(apply(countm1,2,function(x) which.max(x)==2))/(length(Samping[[i]]$test))))
  acc2<-rbind(acc2,c(best,sum(apply(countm2,2,function(x) which.max(x)==6))/(length(Samping[[i]]$test))))
  acc3<-rbind(acc3,c(best,sum(apply(countm3,2,function(x) which.max(x)==10))/(length(Samping[[i]]$test))))
  Best<-c(Best,best)
  cc1<-apply(test[,grep(".6P|X6.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
  cc2<-apply(test[,grep(".7P|X7.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
  cc3<-apply(test[,grep("NC.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
  ccc1=sum(apply(cc1,2,function(x) which.max(x)==2))/5
  ccc2=sum(apply(cc2,2,function(x) which.max(x)==6))/5
  ccc3=sum(apply(cc3,2,function(x) which.max(x)==10))/5
  ccc1
  ccc2
  ccc3
  rlt<-c(mean(Best),ccc1,ccc2,ccc3)
  Rlt<-rbind(Rlt,rlt)
  }
  return(Rlt)
}

######################################################################
##################### One-time test ##################################
######################################################################
# The most easy situation is how the performance for one-time prediction
# Take different `t` and then use all the samples, give the prediction arracy
data<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.txt",head=T,row.names=1,sep="\t")
save(data,file="MHL4.RData")
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
# Response: Random select 5 samples for each group, leave them alone and then apply our modle to these samples.
out1<-grep(".6P|X6.P",colnames(data))[sample(1:30,5)]
out2<-grep(".7P|X7.P",colnames(data))[sample(1:30,5)]
out3<-grep("NC.P",colnames(data))[sample(1:75,5)]
input<-data[match(rownames(bio),rownames(data)),]
acc1<-c()
acc2<-c()
acc3<-c()
Lnum1<-c()
Lnum2<-c()
Lnum3<-c()
Num<-c()
for(j in seq(0,1,0.005)){
counts1<-apply(input[,grep(".6P|X6.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
counts2<-apply(input[,grep(".7P|X7.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
counts3<-apply(input[,grep("NC.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
num<-data.frame(id=j,c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                    c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                    c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
Num<-rbind(Num,num)
}
acc<-sweep(Num, 2,c(1,30,29,75) , `/`)

pdf("one-time-prediction-test-acc.pdf")
plot(c1~id,acc,type="l",col=2,cex=2,lwd=2,ylim=c(0,1),xlab="MHL threshold (t)",ylab="Accuracy of the prediction")
lines(c2~id,acc,lty=2,col=3,cex=2,lwd=2)
lines(c3~id,acc,lty=3,col=4,cex=2,lwd=2)
legend("topright",legend=c("CCP","LCP","NCP"),lty=c(1,2,3),col=c(2,3,4),lwd=2, bty="n")
dev.off()

