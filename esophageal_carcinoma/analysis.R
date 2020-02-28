
setwd("/media/Home_Raid1/shg047/NAS3/Minghua2016/bedgraph")
data<-read.table("data.txt",head=T,row.names=1,sep="\t")
pdf("CG.pdf")
par(las=2,cex.lab=0.5,cex.axis=0.5, cex.lab=0.35,cex.axis=0.35)
boxplot(data,na.rm=T,outline=F,horizontal=T)
dev.off()
Mean<-apply(data,2,function(x) mean(x,na.rm=T))


xx<-data[,match(c("N13.bedgraph","T13.bedgraph"),colnames(data))]

xx<-data[,match(c("N22.bedgraph","N10.bedgraph"),colnames(data))]

chr1:111217330-111217331


pdf("CG.pdf")
boxplot(data,na.rm=T,outline=F)
dev.off()


data<-read.table("CHH.txt",head=T,row.names=1,sep="\t")
pdf("CHH.pdf")
boxplot(data,na.rm=T,outline=F,ylim=c(0,1))
dev.off()

ci95<-function(x){
  error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
  m<-round(mean(x),2)
  d<-round(mean(x)-error,2)
  u<-round(mean(x)+error,2)
  paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
}

data<-read.table("CG.txt",head=T,row.names=1,sep="\t")

pdf("CG.pdf")
par(las=2,cex.lab=0.5,cex.axis=0.5)
boxplot(data,na.rm=T,outline=F,ylim=c(0,1),horizontal=T)
dev.off()

lapply(rownames(data),function(x) substr(x,1,9))

# data input
data<-read.table("data.txt",head=T,row.names=1,sep="\t")
x1<-grep("^N[0-9]",colnames(data))
x2<-grep("^T[0-9]",colnames(data))
colnames(data)[x1]
colnames(data)[x2]

# remove suspected samples
Mean1<-apply(data[,x1],2,function(x) mean(x,na.rm=T))
Mean2<-apply(data[,x2],2,function(x) mean(x,na.rm=T))
Mean<-data.frame(Mean1,Mean2)
which(abs(Mean1-Mean2)<1)
Mean[which(abs(Mean1-Mean2)<1),]

data<-data[,c(x1[which(abs(Mean1-Mean2)>=1)],x2[which(abs(Mean1-Mean2)>=1)])]

#re-analysis
x1<-grep("^N[0-9]",colnames(data))
x2<-grep("^T[0-9]",colnames(data))

Cor<-c()
for(i in 1:length(x1)){
  Cor<-c(Cor,cor(data[,x1[i]],data[,x2[i]],use="pairwise.complete.obs"))  
}

p<-ttestFunction(data,x1,x2)

for(i in 1:50){
t.test(data[i,x1],data[i,x2])$p.value
}

write.table(p,file="Pvalue.txt",col.names=NA,row.names=T,sep="\t",quote=F)

ttestFunction<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],5)
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),sum(! is.na(data[i,x1]))<2,sum(! is.na(data[i,x2]))<2,all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]), na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-as.numeric((mean(as.numeric(data[i,x1]),na.rm=T)-mean(as.numeric(data[i,x2]),na.rm=T)))
      output[i,3]<-"t-test"
      output[i,4]<-mean(as.numeric(data[i,x1]),na.rm=T)
      output[i,5]<-mean(as.numeric(data[i,x2]),na.rm=T)
      print(i)
    }
  }
  out<-cbind(rownames(data),output)
  P.Adj<-p.adjust(output[,1],method="fdr")
  out<-data.frame(out[,1],out[,2],P.Adj,out[,3:6])
  colnames(out)=c("RowName","P-value","P.FDR","Delta","Test","M(G1)","M(G2)")
  return(out)
}