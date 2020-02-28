

setwd("/home/sguo/methylation")
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

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

library("impute")
map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")
pdf("methylation_block.pdf")
par(mfrow=c(3,3))
file=list.files(pattern="*pair.RData")
for(j in 1:length(file)){
  load(file[j])
  cancer<-substr(file[j],1,4)
  data<-RawNARemove(data)
  data<-impute.knn(data)$data
  map<-map[map[,4] %in% rownames(data),]
  newdata<-data[match(map[,4],rownames(data)),]
  if(nrow(map)==nrow(newdata)){
    a<-map[,2]
    i=1
    tmp<-c()
    rlt<-c()
    rowname<-c()
    index<-0
    while(i < length(a)){
      start=a[i]
      end=a[i+1]
      end-start
      if(end-start<100 && end-start>0){
        tmp<-c(tmp,i)
        i=i+1
      }else{
        if(length(tmp)>4){
          index=index+1
          tmp<-c(min(tmp),max(tmp),length(tmp),a[max(tmp)]-a[min(tmp)],round(length(tmp)/(a[max(tmp)]-a[min(tmp)]),4))
          rlt<-rbind(rlt,tmp)
          tmp2<-paste(map[tmp[1],1],":",map[tmp[1],2],"-",map[tmp[2],2],sep="")
          rowname<-c(rowname,tmp2)
        }
        tmp<-c()
        i=i+1
      }
    }
    rownames(rlt)<-rowname
    save(rlt,file=paste(cancer,"mh.rlt.RData",sep="."))
    
    cor<-c()
    for(j in 1:nrow(rlt)){
      cor1<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],seq(1,ncol(newdata),by=2)]),use="complete.obs")) # cancer
      cor2<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],seq(2,ncol(newdata),by=2)]),use="complete.obs")) # normal
      tmp<-c(cor1,cor2)
      cor<-rbind(cor,tmp)
    }
    rownames(cor)<-rowname
    save(cor,file=paste(cancer,"mh.cor.RData",sep="."))
    
    plot(density(cor[,1]),xlim=c(0,1),ylim=c(0,2.5),lwd=3,col="red",main=cancer)
    lines(density(cor[,2]),xlim=c(0,1),ylim=c(0,2.5),lwd=3,col="blue")
    legend("topright",legend=c("cancer","normal"),col=c("red","blue"),lwd=3,lty=1,bty="n")
  }
}


setwd("/home/sguo/methylation")
file=list.files(pattern="*mh.cor.RData")
for(i in 1:length(file)){
  load(file[i])
  cancer=substr(file[i],1,4)
  UCEC.mh.cor.RData
  rownames(cor[which(cor[,2]>0.6),])

}






