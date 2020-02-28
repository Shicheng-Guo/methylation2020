# Takai and Jones's sliding-window algorithm

# chr1  200
# chr1  249
# chr1  290
# chr1  305

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

Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}



map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")



setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod")
load("COADmeth450.pair.RData")
map<-read.table("GPL13534.sort.bed",as.is=T,sep="\t")
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


methcor<-function(data){
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
return(cor)
}

load("COADmeth450.pair.RData")
coadcor<-methcor(data)
load("PRADmeth450.pair.RData")
pradcor<-methcor(data)
load("LUSCmeth450.pair.RData")
lusccor<-methcor(data)


par(mfrow=c(2,2))
plot(density(rlt[,3]),xlim=c(4,15),lwd=3,col="blue",main="")
plot(density(1/rlt[,5]),xlim=c(0,70),lwd=3,col="blue",main="")
plot(density(cor[,1]),xlim=c(0,1),ylim=c(0,2.5),lwd=3,col="red",main="")
lines(density(cor[,2]),xlim=c(0,1),ylim=c(0,2.5),lwd=3,col="blue",main="")
legend("topright",legend=c("cancer","normal"),col=c("red","blue"),lwd=3,lty=1,bty="n")


load("PRADmeth450.pair.RData")

newdata2<-data[match(map[,4],rownames(data)),]
library("impute")
newdata2<-impute.knn(newdata2)$data
cor<-c()
for(j in 1:nrow(rlt)){
  cor1<-mean(cor(t(newdata2[rlt[j,1]:rlt[j,2],seq(1,ncol(newdata),by=2)]),use="complete.obs")) # cancer
  cor2<-mean(cor(t(newdata2[rlt[j,1]:rlt[j,2],seq(2,ncol(newdata),by=2)]),use="complete.obs")) # normal
  tmp<-c(cor1,cor2)
  cor<-rbind(cor,tmp)
}

par(mfrow=c(2,2))
plot(density(rlt[,3]),xlim=c(3,15),lwd=3,col="blue",main="")
plot(density(1/rlt[,5]),xlim=c(0,70),lwd=3,col="blue",main="")

cor=coadcor
plot(density(cor[,1]),xlim=c(0,1),ylim=c(0,2.5),lwd=3,col="red",main="coad")
lines(density(cor[,2]),xlim=c(0,1),ylim=c(0,2.5),lwd=3,col="blue",main="")
legend("topright",legend=c("cancer","normal"),col=c("red","blue"),lwd=3,lty=1,bty="n")

cor=lusccor
plot(density(cor[,1]),xlim=c(0,1),ylim=c(0,3),lwd=3,col="red",main="lusc")
lines(density(cor[,2]),xlim=c(0,1),ylim=c(0,3),lwd=3,col="blue",main="")
legend("topright",legend=c("cancer","normal"),col=c("red","blue"),lwd=3,lty=1,bty="n")


data[order(lusccor,decreasing=T)[1:10],]

cor=pradcor
plot(density(cor[,1]),xlim=c(0,1),ylim=c(0,3.5),lwd=3,col="red",main="")
lines(density(cor[,2]),xlim=c(0,1),ylim=c(0,3.5),lwd=3,col="blue",main="")
legend("topright",legend=c("cancer","normal"),col=c("red","blue"),lwd=3,lty=1,bty="n")




library("GEOquery")
GSE56044 <- getGEO("GSE42861",destdir="/home/sguo/monod/data/geo")
save(GSE56044, file="GSE42861_matrix.Rdata")





