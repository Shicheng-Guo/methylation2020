# This software is Copyright © 2017 The Regents of the University of California. All Rights Reserved.
#  
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
#  
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# invent@ucsd.edu

# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#  
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS.


#!/usr/bin/R
# Guassian distribution based assignment of tissue-of-origin and cancer status.
# Contact: Shicheng Guo
# Version 1.3
# Update: Dec/29/2016
# Input: the ts-MHL counts for each samples in each reference given specific MHL positive threshold

##################################################################################################
############################################# CCP ################################################
##################################################################################################

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

### Figure 
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
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

par(mfrow=c(2,2))
hist(LCP[6,],breaks=20,xlim=c(0,60),col="darkblue",border="darkblue",ylab="Counts of samples",xlab="Counts of ls-MHL",main="Colon")
hist(LCP[11,],breaks=8,xlim=c(0,60),col="darkblue",border="darkblue",ylab="Counts of samples",xlab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)

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

# bar plot show Z-score for one lung sample

barplot(Zmax[,2],col="darkblue",ylab="Z-score")

file=list.file("*txt")
for(i in file){
  d<-read.table(i)
}

# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==6))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==6)))/ncol(Zmax)


p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="bonferroni"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)


## Figure C
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

# barplot of AUC for Lung cancer
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Lung+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")


####################################################
################## Normal ##########################
####################################################
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

par(mfrow=c(3,4))
for(i in 1:(nrow(normal))){
hist(normal[i,],col="darkblue",breaks=30,xlim=c(0,50),xlab=rownames(normal)[i],border ="darkblue",main="",ylab="Counts of rsMHL",cex.axis=1.15,cex.lab=1.5)  
abline(v=30,lty=3,col="red",lwd=2)
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


