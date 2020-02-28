

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


map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")
pdf("correlationbetweencancers.pdf")
par(mfrow=c(3,3))

corcal<-function(data){
data<-RawNARemove(data)
head(map)
map<-map[map[,4] %in% rownames(data),]
a<-map[,2]
i=1
tmp<-c()
rlt<-c()
index<-0
while(i < length(a)){
  start=a[i]
  end=a[i+1]
  end-start
  if(end-start<100){
    tmp<-c(tmp,i)
    i=i+1
  }else{
    if(length(tmp)>4){
      index=index+1
      tmp<-c(min(tmp),max(tmp),length(tmp),a[max(tmp)]-a[min(tmp)],round(length(tmp)/(a[max(tmp)]-a[min(tmp)]),4))
      rlt<-rbind(rlt,tmp)
    }
    tmp<-c()
    i=i+1
  }
}

  newdata<-data[match(map[,4],rownames(data)),]
  library("impute")
  newdata<-impute.knn(newdata)$data
  
  cor<-c()
  for(j in 1:nrow(rlt)){
    cor1<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],seq(1,ncol(newdata),by=2)]),use="complete.obs")) # cancer
    cor2<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],seq(2,ncol(newdata),by=2)]),use="complete.obs")) # normal
    tmp<-c(cor1,cor2)
    cor<-rbind(cor,tmp)
  }

  plot(density(cor[,1]),xlim=c(0,1),ylim=c(0,3.5),lwd=3,col="red",main="")
  lines(density(cor[,2]),xlim=c(0,1),ylim=c(0,3.5),lwd=3,col="blue",main="")
  legend("topright",legend=c("cancer","normal"),col=c("red","blue"),lwd=3,lty=1,bty="n")
}


file=list.files(pattern="*.pair.RData")
for(i in 1:length(file)){
load(file[i])
corcal(data)
print(file[i])
}
dev.off()


