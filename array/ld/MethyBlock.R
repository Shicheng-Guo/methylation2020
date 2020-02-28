setwd("/home/gsc/methylation")

methyblock<-function(){
library("plotrix")
RawNARemove<-function(data,missratio=0.3){
  NaRaw<-which(apply(data,1,function(x) (is.na(x)))==T)
  if(length(NaRaw)>(1-missratio)*dim(data)[2]){
    data1<-data[-NaRaw,]
  }else{
    data1<-data;
  }
  data1
}
ColNARemove<-function(data,missratio=0.3){
  NaCol<-which(apply(data,2,function(x) (is.na(x)))==T)
  if(length(NaCol)>(1-missratio)*dim(data)[1]){
    data1<-data[,-NaCol]
  }else{
    data1<-data;
  }
  data1
}
for (i in 24:1){
  file<-paste("Chr",i,"_450kMerge.txt.trans",sep="")
  data<-read.table(file,head=T,row.names=1,sep="\t",as.is=F)
  data<-t(data)
  data<-RawNARemove(data)
  data<-ColNARemove(data)
  output<-paste("Chr",i,".methy.RData",sep="")
  save(data,file=output)
  cor<-matrix(NA,dim(data)[2],dim(data)[2])
  for(j in 1:dim(data)[2]){
     for(k in j:dim(data)[2]){
  cor[j,k]<-cor(data[,j],data[,k],use="na.or.complete")
  }
  }
  output1<-paste("Chr",i,".cor.RData",sep="")
  save(cor,file=output1)
  output2<-paste("Chr",i,".cor.pdf",sep="")
  print(i)
}
return();
}



