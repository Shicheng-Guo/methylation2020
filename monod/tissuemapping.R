#!/usr/bin/R
# Guassian Distribution based cancer plasma tissue-of-origin and pathology assignment.
# Contact: Shicheng Guo
# Version 1.3
# Update: Dec/29/2016
# Input: the ts-MHL counts for each samples in each reference given specific MHL positive threshold (>0.13)


setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

#################################################
################## CCP ##########################
#################################################

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


p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="fdr"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)
write.table(p1,file="colon.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)
getwd()

# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==2))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))/ncol(Zmax)


### Figure 
par(mfrow=c(4,3),mar=c(4,4,2,1))
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



#################################################
################## LCP ##########################
#################################################

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

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
Zmax

p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="fdr"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)
write.table(p1,file="lung.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)

## Figure C
par(mfrow=c(4,3),mar=c(4,4,2,1))
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




####################################################
################## Normal ##########################
####################################################

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

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
p1<-matrix(p.adjust(p,method="fdr"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)
write.table(p1,file="normal.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)



