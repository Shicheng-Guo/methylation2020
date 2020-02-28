
# Supplementary Figure 11-c
Zmax<-matrix(nrow=nrow(CCP),ncol=ncol(CCP))
for(i in 1:ncol(CCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,60)
    idx2<-idx[which(! idx %in% idx1)]
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    zmp <- (CCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(CCP)
colnames(Zmax)=colnames(CCP)
Zmax
par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i])
}
par(mfrow=c(2,2))
library("beeswarm")
input<-list(CCP=Zmax[2,],WB=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:4,pch=16,method="center",ylab="Z-Score")
p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="bonferroni"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)
sum(p1[2,]>-log(0.05,10))
length(unique(c(which(p1[2,]>-log(0.05,10)),which(p1[11,]>-log(0.05,10)))))
library("gplots")
heatmap.2(p1,Rowv=F,Colv=F,density.info="none",trace="none")
write.table(p1,file="colon.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)
getwd()
# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==2))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))/ncol(Zmax)
par(mfrow=c(4,3),mar=c(4,4,2,1))
AUC<-c()
for(k in 1:11){
  xx<-c()
  sen<-c()
  spe<-c()
  for(fdr in seq(0,1,by=0.01)){
    xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
    sen<-c(sen,xx/30)
    spe<-c(spe,1-fdr)
  }
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main=rownames(Zmax)[k])
  AUC<-c(AUC,mean(sen))
}

Zmax[k,]<-Zmax[2,]+Zmax[11,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Colon+CT")
# barplot of AUC for CRC cancer
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Colon+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")


