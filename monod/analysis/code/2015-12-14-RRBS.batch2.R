#######################################################################################################################
###   Title : Analysis to batch 2 RRBS plasma dataset
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Time :  Sep/23/2015 
###   New: Extract the methylation signals with MethylFreq2Matrix.pl
###   Prerequisite: achieve bed file (target region)
###   Prerequisite: MethylFreq files for all the sample
###   Directory: 512 server: /home/sguo/monod/methyFreq/
#######################################################################################################################

library("impute")
RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
} 

data<-read.table("rrbs.batch2.mhl.txt",row.names=1,head=T,sep="\t",as.is=T,check.names=F)
newdata<-RawNARemove(data)
library("impute")
newdata<-impute.knn(data.matrix(newdata))$data
newdata[1:3,1:3]

newdata <- na.omit(newdata) # listwise deletion of missing
newdata <- scale(newdata) # standardize variables

pdf("hist.plot.pdf")
plot(hist(newdata))
dev.off()




d <- dist(t(newdata), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 

pdf("hclust.pdf")
plot(fit,hang=-1) # display dendogram
dev.off()

library("grDevices")
library("gplots")
pdf("heatmap.batch2.RRBS.pdf")
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(newdata,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.8,labRow=NA)
dev.off()


corsort<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) x))
  bed<-matrix(a,ncol=2,byrow=T)
  bed<-bed[order(bed[,1],as.numeric(bed[,2])),]
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],sep="")})
  return(cor)
}
