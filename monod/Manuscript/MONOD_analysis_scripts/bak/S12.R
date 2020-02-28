#!/usr/bin/R
# Guassian Distribution based cancer plasma tissue-of-origin and pathology assignment.
# Version 1.3
# Update: May/5/2017
# Input: the ts-MHL counts for each samples in each reference given specific MHL positive threshold (>0.13)
####################################################
################## Normal(background) ##########################
####################################################
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))
par(mfrow=c(3,4))
for(i in 1:(nrow(normal))){
  hist(normal[i,],col="darkblue",breaks=30,xlim=c(0,50),xlab=rownames(normal)[i],border ="darkblue",main="",ylab="Counts of rsMHL",cex.axis=1.15,cex.lab=1.5)  
}
Zmax<-matrix(nrow=nrow(normal),ncol=ncol(normal))
for(i in 1:ncol(normal)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,75)
    idx2<-idx[which(! idx %in% idx1)]
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])  
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    
    zmp <- (normal[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(normal)
colnames(Zmax)=colnames(normal)
Zmax
p<-pnorm(-abs(Zmax))
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)

#################################################
################## CCP ##########################
#################################################
par(mfrow=c(2,2))
hist(CCP[2,],breaks=30,xlim=c(0,60),col="darkblue",border="darkblue",xlab="Colon",ylab="Counts of cs-MHL",main="Colon")
hist(CCP[11,],breaks=10,xlim=c(0,60),col="darkblue",border="darkblue",xlab="CT",ylab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)

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
library("beeswarm")
input<-list(Brain=Zmax[1,],Colon=Zmax[2,],Intestine=Zmax[3,],Kidney=Zmax[4,],Liver=Zmax[5,],Lung=Zmax[6,],Pancreas=Zmax[7,],Spleen=Zmax[8,],stomach=Zmax[9,],WBC=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:12,pch=16,method="center",ylab="Z-Score",main="Normal plamsa")

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

#################################################
################## LCP ##########################
#################################################
Zmax<-matrix(nrow=nrow(LCP),ncol=ncol(LCP))
for(i in 1:ncol(LCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,29)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    
    zmp <- (LCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(LCP)
colnames(Zmax)=colnames(LCP)

library("beeswarm")
input<-list(Brain=Zmax[1,],Colon=Zmax[2,],Intestine=Zmax[3,],Kidney=Zmax[4,],Liver=Zmax[5,],Lung=Zmax[6,],Pancreas=Zmax[7,],Spleen=Zmax[8,],stomach=Zmax[9,],WBC=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:12,pch=16,method="center",ylab="Z-Score",main="Normal plamsa")

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
  AUC<-c(AUC,mean(sen))
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main=rownames(Zmax)[k])
}
Zmax[k,]<-Zmax[6,]+Zmax[11,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Lung+CT")

