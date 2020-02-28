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
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("minfiData")
BiocManager::install("missMethyl")
BiocManager::install("minfiData")
BiocManager::install("Gviz")
BiocManager::install("DMRcate")
BiocManager::install("ChAMP")
library("ggplot2")
require("minfi")
library("knitr")
library("limma")
library("minfi")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kmanifest")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("IlluminaHumanMethylationEPICmanifest")
library("RColorBrewer")
library("missMethyl")
library("minfiData")
library("Gviz")
library("DMRcate")
library("stringr")
library("minfi")
library("minfiData")
library("sva")
library("ChAMP")

source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")


RINfun=function(yorig){
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
############################################################################################
### Section 2. read the idat
############################################################################################
baseDir="C:\\Users\\shg047\\Documents\\GitHub\\AtrialFibrillation\\idat"
baseDir="/mnt/bigdata/Genetic/Projects/shg047/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"

setwd(baseDir)
list.files()
dataDirectory <- baseDir
list.files(dataDirectory, recursive = TRUE)
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(base = baseDir, targets = targets)
phenoData <- pData(RGset)
manifest <- getManifest(RGset)
head(getProbeInfo(manifest))
MSet <- preprocessRaw(RGset) 
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
bVals <- getBeta(RSet)
mVals <- getM(RSet)
CN <- getCN(RSet)
sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)
head(mVals[,1:3])
head(bVals[,1:3])
GRset <- mapToGenome(RGset)
annotation <- getAnnotation(GRset)
names(annotation)

qc <- getQC(MSet)
pdf("Figure_S10.pdf")
plotQC(qc)
dev.off()

pdf("Figure_S11.pdf")
densityPlot(MSet, sampGroups = phenoData$Case_Control)
dev.off()

pdf("Figure_S12.pdf")
densityBeanPlot(MSet, sampGroups = phenoData$Case_Control)
dev.off()

pdf("Figure_S13.pdf")
controlStripPlot(RGset, controls="BISULFITE CONVERSION II")
dev.off()

qcReport(RGset, pdf= "Figure_S14_qcReport.pdf")


predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
pdf("Figure_S15_genderprediction.pdf")
plotSex(getSex(GRset, cutoff = -2))
dev.off()

MSet.illumina <- preprocessIllumina(RGset, bg.correct = TRUE,normalize = "controls")
MSet.swan <- preprocessSWAN(RGset)
GRset.quantile <- preprocessQuantile(RGset, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)


MSet.noob <- preprocessNoob(RGset)
GRset.funnorm <- preprocessFunnorm(RGset)
snps <- getSnpInfo(GRset)
head(snps,10)
GRset <- addSnpInfo(GRset)
GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)

library(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(RGset)

beta <- getBeta(GRset.funnorm)
phen  <- pData(GRset.funnorm)$Case_Control
dmp <- dmpFinder(beta, pheno = phen,type="categorical",qCutoff=0.05)
head(dmp)
write.table(dmp,file="AtrialFibrillation.DMR.2.txt",sep="\t",col.names = NA,row.names = T,quote=F)

pheno <- pData(GRset.funnorm)$status
designMatrix <- model.matrix(~ pheno)


dmrs <- bumphunter(GRset.funnorm, design = designMatrix,cutoff = 0.2, B=0, type="Beta")
dmrs <- bumphunter(GRset.funnorm, design = designMatrix,cutoff = 0.2, B=1000, type="Beta")









pdf("Figure_S1.pdf")
densityPlot(RGset,xlim=c(0,1),sampGroups = RGset$Sample_Group,main = "Beta", xlab = "Beta",cex=0.1)
detP <- detectionP(RGset)
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP), col=pal[factor(targets$Case_Control)], las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Case_Control)), fill=pal,bg="white")
dev.off()

qcReport(RGset, sampNames=targets$ID, sampGroups=targets$Case_Control,pdf="Figure_S2.qcReport.pdf")

head(targets)
keep <- colMeans(detP) < 0.05
rgSet <- RGset[,keep]
rgSet
targets <- targets[keep,]
targets[,1:5]
mSetSq <- preprocessQuantile(rgSet) 
mSetRaw <- preprocessRaw(rgSet)

pdf("Figure_S4.pdf")
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Case_Control,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Case_Control)), text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Case_Control,main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Case_Control)), text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("Figure_S5.pdf")
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)],pch=16,cex=1.5)
legend("top", legend=levels(factor(targets$Case_Control)), text.col=pal,bg="white", cex=0.7,pch=16,col=pal)
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Young_Old)],pch=16,cex=1.5)
legend("top", legend=levels(factor(targets$Young_Old)), text.col=pal,bg="white", cex=0.7,pch=16,col=pal)
dev.off()

pdf("Figure_S6.pdf")
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)], dim=c(1,3),pch=16,cex=1.5)
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,cex=0.7, bg="white",pch=16,col=pal)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Case_Control)], dim=c(2,3),pch=16,cex=1.5)
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,cex=0.7, bg="white",pch=16,col=pal)

plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)], dim=c(3,4),pch=16,cex=1.5)
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,cex=0.7, bg="white",pch=16,col=pal)
dev.off()

mSetSqFlt<-mSetSq
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

head(bVals[,1:5])
head(mVals[,1:5])

pdf("Figure_S6.pdf")
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Case_Control, main="Beta values", legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Case_Control)),text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Case_Control, main="M-values",legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Case_Control)), text.col=brewer.pal(8,"Dark2"))
dev.off()

ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
####################################################################################################################################
### Section 3. DMS and DMR analysis
####################################################################################################################################
MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE,normalize = "controls", reference = 2)
mdsPlot(MSet.norm, numPositions = 1000, sampGroups = MSet.norm$Sample_Group, sampNames = MSet.norm$Sample_Name)
mset <- MSet.norm
M <- getM(mset, type = "beta", betaThreshold = 0.001)
dmp <- dmpFinder(M, pheno=mset$Case_Control, type="categorical",qCutoff=0.05)

write.table(dmp,file="AtrialFibrillation.DMR.txt",sep="\t",col.names = NA,row.names = T,quote=F)

pdf("Figure_S6_dmp.pdf")
plotCpg(mset, cpg=rownames(dmp)[1], pheno=mset$Case_Control)
dev.off()
mset <- mset[rownames(dmp),]
mse <- mapToGenome(mset)
rowData(mse)
mcols(rowData(mse)) <- cbind(mcols(rowData(mse)), dmp)

mSetSq <- preprocessQuantile(rgSet) 
mSetRaw <- preprocessRaw(rgSet)

pdf("Figure_S7.pdf")
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)])
legend("top",legend=levels(factor(targets$Case_Control)),text.col=pal, bg="white", cex=0.7)
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Case_Control)])
legend("top",legend=levels(factor(targets$Case_Control)),text.col=pal, bg="white", cex=0.7)
dev.off()

pdf("Figure_S8.pdf")
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Case_Control)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Case_Control)), text.col=pal, cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=pal[factor(targets$Case_Control)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Case_Control)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common",col=pal[factor(targets$Case_Control)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Case_Control)), text.col=pal, cex=0.7, bg="white")
dev.off()



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


