#####################################################################
###   Title : Genome-wide DNA methylation analysis for RA patient
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Section 1. Function predifinition 
###   Section 2. Data Cleaning
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
###   Section 5. GEO Validation (GSE34639,GSE27895)
#####################################################################
####################################################################################################################################
### Section 1. function predifinition 
####################################################################################################################################
library("ggplot2")

RINfun=function(yorig)
{
  yranks=rank(yorig)
  tempp=(yranks-.5)/(length(yranks))
  return(qnorm(tempp))
}
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
####################################################################################################################################
### Section 2. Data Cleaning
####################################################################################################################################
setwd("/home/sguo/Dropbox/Project/methylation/RA/analysis")
sam1<-read.table("sam1.txt",head=T,sep="\t",as.is=T)
sam2<-read.table("sam2.txt",head=T,sep="\t",as.is=T)
sam<-merge(sam2,sam1,by.x=1,by.y=4)
write.table(sam,file="sam.txt",col.names=T,row.names=F,sep="\t",quote=F)

setwd("/home/shg047/meth450/dyh")
data1<-read.table("d0413.txt",head=T,row.names=1,sep="\t",as.is=T)
data2<-read.table("d0514.txt",head=T,row.names=1,sep="\t",as.is=T)
data<-cbind(data1,data2)
sum(is.na(data))
write.table(data,file="rawdata.txt",col.names=NA,row.names=T,sep="\t",quote=F)
####################################################################################################################################
### Section 2. PCA Analysis
####################################################################################################################################
setwd("/home/shg047/meth450/dyh/analysis")
data<-read.table("rawdata.txt",head=T,sep="\t",as.is=T,row.names=1)
head(data)
sam<-read.table("sam.txt",head=T,sep="\t",as.is=T)
head(sam)
data<-data[,match(paste(sam$Sample.Group,"AVG_Beta",sep="."),colnames(data))]
head(sam)
library("impute")
data<-data.matrix(data)
data<-RawNARemove(data)

d<-impute.knn(data)$data
colnames(d)<-colnames(data)
rownames(d)<-rownames(data)

data<-t(d)                                       # requre row is sample and column is variable                           
rownames(data)
phen<-c(rep(1,12),rep(0,12))                     # Ecential
dim(data)

pca <- prcomp(data)
save(pca,file="pca.result.RData")
scores <- data.frame(phen, pca$x[,1:4])

sam[match(sapply(strsplit(rownames(data),"[.]"),function(x) unlist(x)[1]),sam$Sample.Group),]


head(sam)


pdf("component.plot.1.2.pdf")
pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=factor(phen)) + theme(legend.position="none")
pc1.2
dev.off()
pdf("component.plot.1.3.pdf")
pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=factor(phen)) + theme(legend.position="none")
pc1.3
dev.off()
pdf("component.plot.1.4.pdf")
pc1.4 <- qplot(x=PC1, y=PC4, data=scores, colour=factor(phen)) + theme(legend.position="none")
pc1.4
dev.off()
pdf("component.plot.2.3.pdf")
pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=factor(phen)) + theme(legend.position="none")
pc2.3
dev.off()

####################################################################################################################################
### Section 3. Differential Methylation Regions Identification
####################################################################################################################################
setwd("/home/shg047/meth450/dyh/analysis")
data<-read.table("rawdata.txt",head=T,sep="\t",as.is=T,row.names=1)
sam<-read.table("sam.txt",head=T,sep="\t",as.is=T)
data<-data[,match(paste(sam$variable,"AVG_Beta",sep="."),colnames(data))]
# differential methylation test1  # biocLite("outliers")
library("impute")
library("outliers")
data<-data.matrix(data)
data<-RawNARemove(data)
d<-impute.knn(data)$data
colnames(d)<-colnames(data)
rownames(d)<-rownames(data)
data<-d
pvalue<-c()

for(i in 1:nrow(data)){
  x<-data[i,1:12]
  y<-data[i,13:24]
  for(j in 1:3){
  x<-rm.outlier(x, fill = T, median = T, opposite = FALSE)
  y<-rm.outlier(y, fill = T, median = T, opposite = FALSE)
  }
  tmp<-t.test(x,y,pair=T,na.rm=T)
  p1<-tmp$p.value
  t1<-round(tmp$estimate,3)
  u<-mean(x,na.rm=T)
  v<-mean(y,na.rm=T)
  r1<-round(u/v,3)
  tmp1<-c(rownames(data)[i],p1,t1,r1,u,v)
  pvalue<-rbind(pvalue,tmp1)
  print(i)
}
colnames(pvalue)=c("cpg","pvalue","delta","Ratio","RA","HP")
write.table(pvalue,file="differential.analysis.outlier2.txt",col.names=T,row.names=F,sep="\t",quote=F)

pdf("distribution.pdf")
plot(density(as.numeric(data[1,1:12])),col="red")
lines(density(as.numeric(data[1,13:24])),col="blue")
dev.off()
library("Deducer")
perm.t.test(x,y,statistic=c("t","mean"),alternative=c("two.sided", "less", "greater"), midp=TRUE, B=1000000)
x<-rnorm(50,1,1)
y<-rnorm(50,0.8,1)
sum(pvalue<0.05/nrow(data))
data[data<=0.3]<-0
data[data>0.3]<-1
rlt<-apply(data,1,function(x) na.omit(as.numeric(x[1:12]))-na.omit(as.numeric(x[13:24])))
pdf("methylation.freq.change.pdf")
plot(hist(rlt),col="red",main="",xlab="",lwd=2,cex.lab=1.5,xlim=c(-10,10))
dev.off()
####################################################################################################################################
### Section 4. heatmap plot
####################################################################################################################################
setwd("/home/shg047/meth450/dyh/analysis")
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\RA\\analysis")
dms<-read.table("differential.analysis.txt",sep="\t",head=T)
sum(dms$pvalue<0.05/nrow(dms)) # 1202 signficant differential methylation site
data<-read.table("rawdata.txt",head=T,sep="\t",as.is=T,row.names=1)
sam<-read.table("sam.txt",head=T,sep="\t",as.is=T)
data<-data[,match(paste(sam$Sample.Group,"AVG_Beta",sep="."),colnames(data))]
library("heatmap.plus")
library("gplots")
pdf("heatmap.pdf")
heatmapdata<-data[match(dms[which(abs(dms$RA-dms$HP)>0.1 & dms$pvalue<0.05/length(dms$pvalue)),1],rownames(data)),]
colnames(heatmapdata)<-sam[match(sapply(strsplit(colnames(data),"[.]"),function(x) unlist(x)[1]),sam$Sample.Group),1]
heatmap.2(data.matrix(heatmapdata),col=redgreen(7),trace="none",density.info="none")
dev.off()
####################################################################################################################################
### Section 4. GEO Dataset Validation
####################################################################################################################################
setwd("/home/sguo/Dropbox/Project/methylation/RA")

library("GEOquery")
library("affy")
library("simpleaffy")
PhenRecode<-function(dat,p,phen="source_name_ch1"){
  nphen<-match(phen,colnames(p))
  phen<-as.numeric(p[match(colnames(dat),p$geo_accession),nphen])
  phen
}
PCASampleStructure<-function(data=t(dat_GSE27895),phen=phen_GSE27895,filename="GSE27895"){
  rlt<-list()
  pca <- prcomp(data)
  rlt$pca<-pca
  output0<-paste(filename,"pca.result.RData",sep=".")
  output1<-paste(filename,"component.plot.1.2.pdf",sep=".")
  output2<-paste(filename,"component.plot.1.3.pdf",sep=".")
  output3<-paste(filename,"component.plot.1.4.pdf",sep=".")
  output4<-paste(filename,"component.plot.2.3.pdf",sep=".")
  output5<-paste(filename,"component.plot.2.4.pdf",sep=".")
  output6<-paste(filename,"component.plot.3.4.pdf",sep=".")
  save(pca,file=output0)
  scores <- data.frame(phen, pca$x[,1:4])
  pdf(output1)
  plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(phen))[i],col=as.numeric(as.factor(phen))[i],cex=0.8)
  }
  legend("topright",pch=unique(as.numeric(as.factor(phen))),col=unique(as.numeric(as.factor(phen))),legend=unique(phen), cex=0.8)
  dev.off()
}
  
GSE27895 <- getGEO("GSE27895",destdir="/home/shg047/meth450/dyh/analysis")
save(GSE27895, file="GSE27895_matrix.Rdata")

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\RA\\analysis")
setwd("/home/sguo/Dropbox/Project/methylation/RA")
load("GSE27895_matrix.Rdata")
dat_GSE27895 <- as.data.frame(exprs(GSE27895[[1]]))
p_GSE27895 <- pData(phenoData(GSE27895[[1]]))
table(p_GSE27895$characteristics_ch1.1) # GSE42752:  11 RA and 12 Normal CD4+
phen_GSE27895<-PhenRecode(dat_GSE27895,p_GSE27895,phen="characteristics_ch1.1")
source("BMIQ_1.3.R")

dim(dat_GSE27895)

BMIQ(beta.v,design.v,nL=3,doH=TRUE,nfit=5000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1)
  

hit<-apply(dat_GSE27895,1,function(x) t.test(x[1:11],x[12:23])$p.value)
for(i in seq(1000,7000,by=1000)){
data<-t(dat_GSE27895[order(hit,decreasing=F)[1:i],])
file=paste("GSE27895",i,sep=".")
rlt<-PCASampleStructure(data=data,phen=phen_GSE27895,filename=file)
}

GSE46650 <- getGEO("GSE46650",destdir="/home/sguo/Dropbox/Project/methylation/RA/analysis")
GSE46650 <- getGEO("GSE46650",destdir="/home/shg047/meth450/dyh/analysis")
save(GSE46650, file="GSE46650_matrix.Rdata")
load("GSE46650_matrix.Rdata")
dat_GSE46650 <- as.data.frame(exprs(GSE46650[[1]]))
p_GSE46650 <- pData(phenoData(GSE46650[[1]]))
table(p_GSE27895$characteristics_ch1.1) # GSE42752:  11 RA and 12 Normal CD4+
phen_GSE27895<-PhenRecode(dat_GSE46650,p_GSE46650,phen="characteristics_ch1.1")



data<-t(dat_GSE27895)     
pca <- prcomp(data,center=T,scale = F)  # Here, input file: row is individual and column is variable
pdf("Figure2.PCA.loading.pdf")
plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
dev.off()
var<-c()
for(i in 1:length(pca$sdev)){
  var[i]<-sum((pca$sdev[1:i])^2)/sum((pca$sdev)^2)
}
pdf("Figure2.PCA.loading.total.increasing.pdf")
plot(var,ylab="total variance",xlab="number of principle components",lwd=2)
dev.off()
save(pca,file="pca.result.RData")

scores <- data.frame(phen_GSE27895, pca$x[,1:3])
col = as.numeric(as.factor(phen_GSE27895))
pheno=phen_GSE27895
pdf("Figure2.PC12.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
for(i in 1:length(scores$PC1)){
  points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.8,lwd=2)
}
dev.off()
pdf("Figure2.PC13.pdf")
plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
for(i in 1:length(scores$PC1)){
  points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
}
dev.off()
pdf("Figure2.PC23.pdf")
plot(x=scores$PC2,y=scores$PC3, xlim=c(min(scores$PC2),max(scores$PC2)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC2",ylab="PC3")
for(i in 1:length(scores$PC1)){
  points(scores$PC2[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
}
dev.off()



