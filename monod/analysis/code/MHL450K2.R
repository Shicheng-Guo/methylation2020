

# Takai and Jones's sliding-window algorithm
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

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) strsplit(x,"-")))
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
dev.off()



par(mfrow=c(2,2))
plot(density(rlt[,3]),xlim=c(3,15),lwd=3,col="blue",main="")
plot(density(1/rlt[,5]),xlim=c(0,70),lwd=3,col="blue",main="")


setwd("/home/sguo/methylation")
file=list.files(pattern="*mh.cor.RData")
pdf("cancer-normal.correlation.pdf")
par(mfrow=c(3,3))
for(i in 1:length(file)){
load(file[i])
cancer=substr(file[i],1,4)
plot(cor[,2],cor[,1],cex=0.4,col="blue",xlim=c(0.15,1),ylim=c(0.15,1),xlab="Normal",ylab="Cancer",cex.lab=1.5,cex.axis=1.5,main=cancer)
lines(seq(0,1,by=0.01),seq(0,1,by=0.01),lty=1,lwd=2.5,col="red")
}
dev.off()


setwd("/home/sguo/methylation")
file=list.files(pattern="*mh.rlt.RData")
pdf("high-density.region.pdf")
par(mfrow=c(4,4))
for(i in 1:length(file)){
  load(file[i])
  cancer=substr(file[i],1,4)
  plot(density(rlt[,3]),xlim=c(3,15),lwd=3,col="blue",main=cancer)
  plot(density(1/rlt[,5]),xlim=c(0,70),lwd=3,col="blue",main=cancer)
}
dev.off()



setwd("/home/sguo/methylation")
file=list.files(pattern="*mh.rlt.RData")
par(mfrow=c(4,4))
for(i in 1:length(file)){
  load(file[i])
  cancer=substr(file[i],1,4)
  # print(quantile(rlt[,4]/rlt[,3],c(0.25,0.5,0.75)))
  print(quantile(rlt[,3]))
  
}


# shared high HDRC regions based on normal samples
setwd("/home/sguo/methylation")
file=list.files(pattern="*mh.cor.RData")
load(file[i])
cancer=substr(file[1],1,4)
edge=0.6
newcor<-cor[which(cor[,2]>edge),]
bed1<-cor2bed(rownames(newcor))
print(nrow(bed1))
number<-c()
region<-c(rownames(newcor))
for(i in 2:length(file)){
  load(file[i])
  cancer=substr(file[i],1,4)
  newcor<-cor[which(cor[,2]>edge),]
  bed2<-cor2bed(rownames(newcor))
  bed1<-Rbedtools(functionstring="intersectBed",bed1,bed2,opt.string="-wa")
  ntmp<-c(nrow(bed1),nrow(bed2))
  number<-c(number,nrow(bed1))
  region<-c(region,rownames(newcor))
}

regionx1<-unique(region)
regionx2<-cor2bed(regionx1)
write.table(regionx1,file="conservative.methylation.block.0.6.normal.cor",sep="\t",col.names=F,row.names=F,quote=F)
write.table(regionx2,file="conservative.methylation.block.0.6.normal.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(bed1,file="conservative.methylation.block.C11.0.6.normal.bed",sep="\t",col.names=F,row.names=F,quote=F)

number4<-c(3734,1926,1679,1220,965,878,795,755,665,654,654)
number5<-c(2260,999, 880, 655, 563, 471, 415, 401, 366, 357, 357)
number6<-c(1181,529, 480, 350, 291, 211, 184, 176, 164, 155, 155)

plot(number4,type="l",lwd=3,col="red",ylim=c(150,3000),xlab="Numbers of cancer",ylab="Number of region",cex.axis=1.25,cex.lab=1.25)
lines(number5,lwd=3,col="blue")
lines(number6,lwd=3,col="black")
legend("topright",legend=c("cut-off=0.6","cut-off=0.5","cut-off=0.4"),col=c("red","blue","black"),lwd=3,lty=1,bty="n",cex=1.25)


region1<-table(region)
region2<-unique(region)

p=length(region)/(length(region2)*11)
pvalue<-c()
for(i in 1:length(table(region))){
ptmp<-binom.test(table(region)[i],11,p,alternative="two.sided")$p.value
pvalue<-c(pvalue,ptmp)
}

sum(pvalue<0.05/length(pvalue))



number4<-rev(c(3734,1926,1679,1220,965,878,795,755,665,654,654))
number5<-rev(c(2260,999, 880, 655, 563, 471, 415, 401, 366, 357, 357))
number6<-rev(c(1181,529, 480, 350, 291, 211, 184, 176, 164, 155, 155))
number4<-rev(cumsum(number4)/sum(number4))
number5<-rev(cumsum(number5)/sum(number5))
number6<-rev(cumsum(number6)/sum(number6))
plot(number4,type="l",lwd=3,lty=2,col="red",ylim=c(0,1),xlab="Numbers of cancer",ylab="proportion of region",cex.axis=1.25,cex.lab=1.25)
lines(number5,lwd=3,lty=2,col="blue")
lines(number6,lwd=3,lty=2,col="black")
legend("topright",legend=c("cut-off=0.6","cut-off=0.5","cut-off=0.4"),col=c("red","blue","black"),lwd=3,lty=3,bty="n",cex=1.25)





# shared high HDRC regions based on cancer samples
setwd("/home/sguo/methylation")
file=list.files(pattern="*mh.cor.RData")
load(file[i])
cancer=substr(file[1],1,4)
edge=0.5
newcor<-cor[which(cor[,2]>edge),]
bed1<-cor2bed(rownames(newcor))
print(nrow(bed1))
number<-c()
region<-c(rownames(newcor))
for(i in 2:length(file)){
  load(file[i])
  cancer=substr(file[i],1,4)
  newcor<-cor[which(cor[,1]>edge),]
  bed2<-cor2bed(rownames(newcor))
  bed1<-Rbedtools(functionstring="intersectBed",bed1,bed2,opt.string="-wa")
  ntmp<-c(nrow(bed1),nrow(bed2))
  number<-c(number,nrow(bed1))
  region<-c(region,rownames(newcor))
}

regionx1<-unique(region)
regionx2<-cor2bed(regionx1)
write.table(regionx1,file="conservative.methylation.block.0.6.cancer.cor",sep="\t",col.names=F,row.names=F,quote=F)
write.table(regionx2,file="conservative.methylation.block.0.6.cancer.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(bed1,file="conservative.methylation.block.C11.0.6.cancer.bed",sep="\t",col.names=F,row.names=F,quote=F)

number41<-(c(2040,1880,1680,1465,1423,1366,1341,1214,989,983))
number51<-(c(1246,1172,1005,925,899,872,864,773,547,536))
number61<-(c(834,792,655,633,603,596,544,397,245,243))
number42<-rev(cumsum(number4)/sum(number4))
number52<-rev(cumsum(number5)/sum(number5))
number62<-rev(cumsum(number6)/sum(number6))

plot(number4,type="l",lwd=3,lty=2,col="red",ylim=c(200,2200),xlab="Numbers of cancer",ylab="proportion of region",cex.axis=1.25,cex.lab=1.25)
lines(number5,lwd=3,lty=2,col="blue")
lines(number6,lwd=3,lty=2,col="black")
legend("topright",legend=c("cut-off=0.6","cut-off=0.5","cut-off=0.4"),col=c("red","blue","black"),lwd=3,lty=3,bty="n",cex=1.25)


plot(number42,type="l",lwd=3,lty=2,col="red",ylim=c(0.2,1),xlab="Numbers of cancer",ylab="proportion of region",cex.axis=1.25,cex.lab=1.25)
lines(number52,lwd=3,lty=2,col="blue")
lines(number62,lwd=3,lty=2,col="black")
legend("topright",legend=c("cut-off=0.6","cut-off=0.5","cut-off=0.4"),col=c("red","blue","black"),lwd=3,lty=3,bty="n",cex=1.25)







