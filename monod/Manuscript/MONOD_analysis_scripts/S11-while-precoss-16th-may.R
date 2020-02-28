#!/usr/bin/R
# Add tumor-specific biomarker to our 10 normal tissue reference and then give distribution.
# Version 1.3
# Update: May/5/2017
# Input: the ts-MHL counts for each samples in each reference given specific MHL positive threshold (>0.13)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("impute")  
#install.packages("gplots")
#install.packages("RColorBrewer")
#install.packages("grDevices")

library("gplots")
library("RColorBrewer")
library("grDevices")
library("impute")

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
}

bed2cor<-function(bed){
  cor<-apply(bed,1,function(x) paste(x[1],":",x[2],"-",x[3],sep=""))
  cor<-gsub(" ","",cor)
  return(cor)
}

gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  GSI<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    if(length(gmax)<1){print(data[i,])}
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    GSI<-c(GSI,gsit)
    #print(paste(gmax,gsit,sep=" "))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=GSI)
  return(rlt)
}

TopGSIByCategory<-function(gsi,top=150){
  GSIRlt<-c()
  group<-names(table(gsi$group))
  rank<-c(rep(top,length(group)))
  for (i in 1:length(group)){
    subset=gsi[which(gsi$group==group[i]),]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
  }
  return(na.omit(GSIRlt))
}
topgsi2bio<-function(topgsi){
  cor2bed<-function(cor){
    cor<-as.character(cor)
    a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
    bed<-matrix(a,ncol=3,byrow=T)
    return(data.frame(bed))
  }
  bio<-data.frame(cor2bed(topgsi[,1]),topgsi[,2:3])
  rownames(bio)<-topgsi[,1]
  return(bio)
}
ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}   
rename<-function(data){
  Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P",colnames(data))]
  colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
  colnames(Data)<-gsub("[.]","-",colnames(Data))
  colnames(Data)[grep("age|new|centenarian",colnames(Data))]<-"WBC"
  colnames(Data)[grep("X7.T|X6.T|SRX|CTT",colnames(Data))]<-"CT"
  colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
  colnames(Data)[grep("X6.P|RRBS.6P",colnames(Data))]<-"CCP"
  colnames(Data)[grep("X7.P|RRBS.7P",colnames(Data))]<-"LCP"
  colnames(Data)[grep("NC.P",colnames(Data))]<-"NCP"
  return(Data)
}

setwd("/home/shg047/oasis/monod/hapinfo")
saminfo<-read.table("/home/shg047/oasis/monod/saminfo/N37Salk.saminfo",sep="\t")
#data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
# feature reduction (WGBS+RRBS missing<60, each plasma category(colon,lung and normal)<60%, Plasma missing<50%)
rm1<-which(apply(data[,grep("X7.P|X6.P|.6P|.7P|NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.5))
rm2<-which(apply(data[,grep("X6.P|.6P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.5))
rm3<-which(apply(data[,grep("X7.P|.7P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.5))
rm4<-which(apply(data[,grep("NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.5))
rm5<-which(apply(data,1,function(x) sum(is.na(x))/length(x)>0.5))
rm<-unique(c(rm1,rm2,rm3,rm4,rm5))
data<-data[-rm,]

#input<-data[match(rownames(bio),rownames(data)),]
# copy our biomarker form supplementary table to: /home/sguo/Dropbox/Project/methylation/monod/Manuscript/MONOD_analysis_scripts/biomarker2.txt (download from supplementary)
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)<-gsub("[.]","-",colnames(Data))
colnames(Data)[grep("age|new|centenarian|WB",colnames(Data))]<-"NT"
colnames(Data)[grep("X7.T|X6.T|SRX|CTT",colnames(Data))]<-"CT"
colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-"NT"
#colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
# for tissue-specific biomarkers

library("impute")
library("preprocessCore")
DATA<-Data[,grep("NT|CT",colnames(Data))] 
#for(i in names(table(colnames(DATA)))){
#  DATA[,colnames(DATA)==i]<-normalize.quantiles(data.matrix(DATA[,colnames(DATA)==i]))
#}
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
DATA<-RawNARemove(DATA,missratio=0.60)
#DATA1<-ColNARemove(DATA,missratio=0.79)
#DATA2<-impute.knn(data.matrix(DATA1))$data
#gsirlt<-gsi(DATA2)
gsirlt<-gsi(DATA)
gsirlt<-subset((gsirlt),group=="CT")
input<-DATA[match(subset((gsirlt),group=="CT")[,1],rownames(DATA)),grep("NT",colnames(DATA))]

rlt<-c()
for(threshold1 in seq(0,0.1,0.01)){
filterout<-rownames(input[which(apply(input,1,function(x) mean(x,na.rm=T)>threshold1)),])
Gsirlt<-gsirlt[-match(filterout,gsirlt[,1]),]
bio2<-topgsi2bio(TopGSIByCategory(Gsirlt,300))
colnames(bio2)<-colnames(bio)
bio3<-rbind(bio,subset((bio2),V5=="CT"))
bio3$V5<-as.character(bio3$V5)
#write.table(data.frame(bio3),"biomarker3.txt",sep="\t",row.names=F,col.names=F,quote=F)
#### Build prediction matrix
#data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
#Data=data[,grep("X7.P|X6.P|.6P|.7P|NC.P",colnames(data))]
Data1<-data[match(rownames(bio3),rownames(data)),]
###########################################################
################## Prediction #############################
###########################################################
# threshold1=0.01 (cancer biomarker threshold in normal tissues)
# threshold1=0.13 (cancer biomarker cancer plasma samples)
CCPF(cc1,threshold1=0.01,threshold2=0.3,zz=1)
LCPF(cc2,threshold1=0.01,threshold2=0.3,zz=1)

zz<-0
for(threshold2 in seq(0,0.2,0.01)){
zz=zz+1
threshold1=0.01       # cancer biomarker threshold in normal tissues
threshold2=0.13       # This threshold can not like the threshold of Figure 5, since there will have remove WBC background
                      # optimize the threshold to achieve the best prediction accuray 
cc<-apply(Data1[,grep("NC.P",colnames(Data1))],2,function(x) tapply(x,bio3$V5,function(x) sum(x>threshold,na.rm=T)))
cc1<-apply(Data1[,grep("X6.P|RRBS.6P",colnames(Data1))],2,function(x) tapply(x,bio3$V5,function(x) sum(x>threshold2,na.rm=T)))
cc2<-apply(Data1[,grep("X7.P|RRBS.7P",colnames(Data1))],2,function(x) tapply(x,bio3$V5,function(x) sum(x>threshold2,na.rm=T)))
CCPF(cc1,threshold1,threshold2,zz)
LCPF(cc2,threshold1,threshold2,zz)
xx1<-mean(cc1[3,]-cc[3,])
xx2<-mean(cc2[3,]-cc[3,])
rlt<-rbind(rlt,c(threshold1,threshold2,xx1,xx2,mean(cc1[11,]),mean(cc2[11,])))
print(c(threshold1,threshold2,xx1,xx2,mean(cc1[11,]),mean(cc2[11,])))
}
}
filename=paste("normal.background",threshold,"jpeg",sep=".")
jpeg(filename)
par(mfrow=c(3,4))
for(i in 1:(nrow(cc))){
  hist(cc[i,],col="darkblue",breaks=30,xlim=c(0,50),xlab=rownames(cc)[i],border ="darkblue",main="",ylab="Counts of rsMHL",cex.axis=1.15,cex.lab=1.5)  
}
dev.off()
write.table(data.frame(cc),"normal-ref.txt",sep="\t",row.names=T,col.names=NA,quote=F)
## For Colon Cancer Plasma
cc1<-apply(Data1[,grep("X6.P|RRBS.6P",colnames(Data1))],2,function(x) tapply(x,bio3$V5,function(x) sum(x>threshold,na.rm=T)))
write.table(data.frame(cc1),"CCP-1.txt",sep="\t",row.names=T,col.names=NA,quote=F)
## For Lung Cancer Plasma
cc2<-apply(Data1[,grep("X7.P|RRBS.7P",colnames(Data1))],2,function(x) tapply(x,bio3$V5,function(x) sum(x>threshold,na.rm=T)))
write.table(data.frame(cc2),"LCP-1.txt",sep="\t",row.names=T,col.names=NA,quote=F)
## Z-Score Matrix
#setwd("C:\\Users\\shg047\\Dropbox\\Project\\methylation\\monod\\Manuscript\\MONOD_analysis_scripts\\bak")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

jpeg("S11-V.jpg")
par(mfrow=c(3,4))
for(i in 1:(nrow(normal))){
  hist(normal[i,],col="darkblue",breaks=30,xlim=c(0,50),xlab=rownames(normal)[i],border ="darkblue",main="",ylab="Counts of rsMHL",cex.axis=1.15,cex.lab=1.5)  
}
dev.off()
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
rownames(p)<-rownames(p)
colnames(p)<-colnames(p)
p1<--log(p,10)

#################################################
################## CCP ##########################
#################################################
library("beeswarm")
CCPF<-function(CCP,threshold1,threshold2,zz){
par(mfrow=c(2,2))
hist(CCP[2,],breaks=30,xlim=c(0,60),col="darkblue",border="darkblue",xlab="Colon",ylab="Counts of cs-MHL",main="Colon")
hist(CCP[3,],breaks=10,xlim=c(0,60),col="darkblue",border="darkblue",xlab="CT",ylab="Counts of cs-MHL",main="Cancer Tissue")
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

write.table(Zmax,file=paste("CCP.ZScore.t1",threshold1,"t2",threshold2,zz,".txt",sep="-"),col.names=NA,row.names=T,sep="\t",quote=F)
jpeg(paste(zz,"ccp.beeswarm.t1",threshold1,"t2",threshold2,".jpg",sep="-"))
rltinput<-list(Brain=Zmax[1,],Colon=Zmax[2,],Intestine=Zmax[4,],Kidney=Zmax[5,],Liver=Zmax[6,],Lung=Zmax[7,],Pancreas=Zmax[8,],Spleen=Zmax[9,],stomach=Zmax[10,],WBC=Zmax[11,],CT=Zmax[3,])
beeswarm(rltinput,col = 2:12,pch=16,method="center",ylab="Z-Score",main="Colon cancer plamsa")
dev.off()
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
Zmax[k,]<-Zmax[2,]+Zmax[3,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Colon+CT")
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Colon+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")
}


#################################################
################## LCP ##########################
#################################################
LCPF<-function(LCP,threshold1,threshold2,z){
  par(mfrow=c(2,2))
  hist(LCP[7,],breaks=30,xlim=c(0,60),col="darkblue",border="darkblue",xlab="Colon",ylab="Counts of cs-MHL",main="Lung")
  hist(LCP[3,],breaks=10,xlim=c(0,60),col="darkblue",border="darkblue",xlab="CT",ylab="Counts of cs-MHL",main="Cancer Tissue")
  Zmax<-matrix(nrow=nrow(LCP),ncol=ncol(LCP))
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
  write.table(Zmax,file=paste("LCP.ZScore.t1",threshold1,"t2",threshold2,zz,".txt",sep="-"),col.names=NA,row.names=T,sep="\t",quote=F)
  jpeg(paste(zz,"Lcp.beeswarm.t1",threshold1,"t2",threshold2,".jpg",sep="-"))
  rltinput<-list(Brain=Zmax[1,],Colon=Zmax[2,],Intestine=Zmax[4,],Kidney=Zmax[5,],Liver=Zmax[6,],Lung=Zmax[7,],Pancreas=Zmax[8,],Spleen=Zmax[9,],stomach=Zmax[10,],WBC=Zmax[11,],CT=Zmax[3,])
  beeswarm(rltinput,col = 2:12,pch=16,method="center",ylab="Z-Score",main="Lung cancer plamsa")
  dev.off()
  
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
  Zmax[k,]<-Zmax[3,]+Zmax[7,]
  xx<-c()
  sen<-c()
  spe<-c()
  for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
  }
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Lung+CT")
}
###########################################################
################## Heatmap #############################
###########################################################
# Take Z-Score Matrix as the Heatmap input matrix
library("gplots")
library("RColorBrewer")
library("grDevices")
library("impute")
bio3<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker3.txt",head=F)  # Download from Supplementary Table
mydata<-data[,grep("X6.P|.6P|X7.P|.7P|NC.P",colnames(data))]
mydata<-mydata[na.omit(match(bed2cor(bio3[,1:3]),rownames(mydata))),]
mydata<-rename(mydata)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")
# perform clustering on rows and columns
cl.row <- hclustfunc(distfunc(mydata))
cl.col <- hclustfunc(distfunc(t(mydata)))
# extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
gr.row <- cutree(cl.row, 6)
# gr.col <- cutree(cl.col, 5)
# require(RColorBrewer)
# col1 <- brewer.pal(6, "Set1")     # the maximum of brewer.pal is 12
bio3<-bio3[na.omit(match(rownames(mydata),bed2cor(bio3[,1:3]))),]
mydata<-mydata[na.omit(match(bed2cor(bio3[,1:3]),rownames(mydata))),]
mydata<-mydata[,order(colnames(mydata))]
colnames(mydata)<-unlist(lapply(colnames(mydata),function(x) unlist(strsplit(x,"[.]"))[1]))
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
#col2<-tol21rainbow[as.numeric(as.factor(cl.col$labels))]
col1 <-tol21rainbow[as.numeric(as.factor(bio3[,4]))]   # the maximum of brewer.pal is 12
col2 <-tol21rainbow[as.numeric(as.factor(colnames(mydata)))]   # the maximum of brewer.pal is 12
col=colorRampPalette(c("yellow", "blue"))(20) 
require(gplots)    
pdf("Figure3A-3-misswithte.pdf")
par(mar=c(5,5,5,0))
heatmaprlt<-heatmap.2(as.matrix(mydata),hclustfun=hclustfunc, distfun=distfunc,
                      RowSideColors=col1, 
                      ColSideColors=col2,
                      labRow=F,
                      trace="none",
                      Colv=F,
                      Rowv=F,
                      col=col,
                      na.color="white",
                      density.info="none")
legend=unique(data.frame(col2,cl.col$labels))
legend(x=0.85,y=0.8,legend=legend[,2],col=as.character(legend[,1]),pch=15,cex = 0.5)
dev.off()


# transfer to genome-miner (make better plot)
setwd("/media/Home_Raid1/shg047/work/monod/hapinfo")
colnames(mydata)<-unlist(lapply(colnames(mydata),function(x) unlist(strsplit(x,"[.]"))[1]))
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
#col2<-tol21rainbow[as.numeric(as.factor(cl.col$labels))]
col1 <-tol21rainbow[as.numeric(as.factor(bio3[,4]))]   # the maximum of brewer.pal is 12
col2 <-tol21rainbow[as.numeric(as.factor(colnames(mydata)))]   # the maximum of brewer.pal is 12
load("mydata.RData")
load("col1.RData")
load("col2.RData")
mydata[mydata>0.13]<-1
mydata[mydata<0.13]<-0
col=colorRampPalette(c("yellow", "blue"))(20) 
require(gplots)    
pdf("Figure3A-3-misswithte-binary2.pdf")
par(mar=c(5,5,5,0))
heatmaprlt<-heatmap.2(as.matrix(mydata),hclustfun=hclustfunc, distfun=distfunc,
                      RowSideColors=col1, 
                      ColSideColors=col2,
                      labRow=F,
                      trace="none",
                      Colv=F,
                      Rowv=F,
                      col=col,
                      na.color="white",
                      density.info="none")
#legend=unique(data.frame(col2,cl.col$labels))
#legend(x=0.85,y=0.8,legend=legend[,2],col=as.character(legend[,1]),pch=15,cex = 0.5)
dev.off()

