setwd("/home/sguo/score");
source("fclust.R")
library("skmeans")
library("plotrix")
library("knitr")
library("fpc")


for (i in 24:1){
file<-paste("Chr",i,"_450kMerge.txt.trans",sep="")
data<-read.table(file,head=T,row.names=1,sep="\t",as.is=F)
output<-paste("Chr",i,".methy.RData",sep="")
save(data,file=output)
cor<-cor(data)
output1<-paste("Chr",i,".cor.RData",sep="")
save(cor,file=output1)
output2<-paste("Chr",i,".cor.pdf",sep="")
pdf(output2)
color2D.matplot(cor,show.values=1,show.legend=T,main="High correlation Region")
dev.off()
print(i)
}


library("plotrix")
library("skmeans")
library("knitr")

x<-matrix(rnorm(1024),nrow=32)
# simulate a correlation matrix with values -0.5 to 0.5
x<-rescale(x,c(-0.5,0.5))
# add a column with the extreme values (-1,1) to calculate
# the colors, then drop the extra column in the result
cellcol<-color.scale(cbind(x,c(-1,rep(1,31))),c(0,1),0,c(1,0))[,1:32]
color2D.matplot(x,cellcolors=cellcol,main="Blue to red correlations")
# do the legend call separately to get the full range
color.legend(0,-4,10,-3,legend=c(-1,-0.5,0,0.5,1),rect.col=color.scale(c(-1,-0.5,0,0.5,1),c(0,1),0,c(1,0)),align="rb")

# Subfunction

ClustValidtion<-function(ClusteResult,ClinicNumericFact){
  rlt<-list()
  data<-cbind(ClinicNumericFact,ClusteResult)
  data2<-RawNARemove(data)
  ClinicNumericFact<-data2[,1:dim(ClinicNumericFact)[2]]
  ClusteResult<-data2[,(dim(ClinicNumericFact)[2]+1):(dim(ClinicNumericFact)[2]+dim(ClusteResult)[2])]
  CorrectRand<-Vi<-matrix(NA,dim(ClusteResult)[2],dim(ClinicNumericFact)[2])
  for (i in 1:dim(ClusteResult)[2]){
    for (j in 1:dim(ClinicNumericFact)[2]){
      result<-cluster.stats(NULL,ClinicNumericFact[,j],ClusteResult[,i],compareonly=T)
      CorrectRand[i,j]<-result$corrected.rand
      Vi[i,j]<-result$vi
    }
  }
  rlt$correctRand<-CorrectRand
  rlt$vi<-Vi
  rlt
}


skmean<-function(data, classnum){
  # suitable to big data especially for col > row matrix, defort to cluste matrix to maximun 10 clusters
  cluster<-matrix(NA, dim(data)[1], classnum-1)
  for (i in 2:classnum){
    ask<-skmeans(data,i)
    file<-paste(data,".sk",i,".RData",sep="")
    cluster[,i-1]<-ask$cluster
    save(ask, file=file)
    print(i)
  }
  cluster
}
RawNARemove<-function(data){
  NaRaw<-which(apply(data,1,function(x) any(is.na(x)))==T)
  if(length(NaRaw)>0){
    data1<-data[-NaRaw,]
  }else{
    data1<-data;
  }
  data1
}

ColNARemove<-function(data){
  NaCol<-which(apply(data,2,function(x) any(is.na(x)))==T)
  if(length(NaCol)>0){
    data1<-data[,-NaCol]
  }else{
    data1<-data;
  }
  data1
}



