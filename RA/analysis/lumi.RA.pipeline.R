#####################################################################
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 7/14/2015
### 
#####################################################################
source("http://bioconductor.org/biocLite.R")
biocLite("lumi")
biocLite("FDb.InfiniumMethylation.hg19")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

setwd("/home/sguo/dyh/idat")
library("ChAMP")
library("lumi")
source("BMIQ_1.3.R")
library("FDb.InfiniumMethylation.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")


# import idat with sample id
sampleInfo<-list.files(pattern="*.idat")
sampleInfo<-gsub("_Red.idat","",sampleInfo)
sampleInfo<-gsub("_Grn.idat","",sampleInfo)
sampleInfo<-unique(sampleInfo)
idat<-importMethyIDAT(sampleInfo, dataPath = getwd(), lib ="FDb.InfiniumMethylation.hg19")
save(idat,file="idat.RData")
# Adjust background level of Illumina Infinium methylation data
idat.bjc <- lumiMethyB(idat, method="bgAdjust2C", separateColor=TRUE)
save(idat.bjc,file="idat.bjc.RData")
# Color bias adjust of Illumina Infinium methylation data
idat.cjc <- lumiMethyC(idat.bjc)
save(idat.cjc,file="idat.cjc.RData")

# Normalize the Illumina Infinium methylation data
idat.qtl <- lumiMethyN(idat.cjc, method="quantile")
save(idat.qtl,file="idat.qtl.RData")

# prepare file upload to GEO with M-value, be sure the first line of GEOSampleInfoTemplate.txt should be start as #
produceGEOSampleInfoTemplate(idat.qtl, fileName = "GEOSampleInfoTemplate.txt") 
produceMethylationGEOSubmissionFile(idat.qtl, idat, sampleInfo="GEOSampleInfoTemplate.txt")

saminfo<-read.table("GEOSampleInfoTemplate.txt",head=T,sep="\t",as.is=T)
write.table(saminfo,file="GEOSampleInfoTemplate.2.txt",sep="\t",quote=F,row.names=F,col.names=T)

load("idat.bjc.RData")
load("idat.cjc.RData")
load("idat.qtl.RData")


# extract raw Beta-value
raw<-matrix(nrow=485577,ncol=48)
dataMatrix1 <- exprs(idat)
beta1<-m2beta(dataMatrix1) # m<-beta2m(beta)
dectP1 <- detection(idat)


# extract M-value and Beta-value
dataMatrix <- exprs(idat.qtl)
beta<-m2beta(dataMatrix) # m<-beta2m(beta)
methy <- assayDataElement(idat.qtl, "methylated")
unmethy <- assayDataElement(idat.qtl, "unmethylated")
dectP <- detection(idat.qtl)

raw<-matrix(nrow=nrow(dataMatrix),ncol=ncol(beta)*2)
ncol<-ncol(beta)
for(i in 1:ncol){
  raw[,2*i-1]=beta1[,i]
  raw[,2*i]=dectP1[,i]
}
colnames(raw)[seq(1,2*ncol,by=2)]<-sam[match(colnames(beta),sam[,2]),1]
colnames(raw)[seq(2,2*ncol,by=2)]<-"Detection Pval"
rownames(raw)<-rownames(beta)

write.table(raw,file="Matrix_processed.txt",sep="\t",quote=F,row.names=T,col.names=NA)

raw<-matrix(nrow=nrow(dataMatrix),ncol=ncol(beta)*3)
ncol<-ncol(beta)
for(i in 1:ncol){
  raw[,3*i-2]=unmethy[,i]
  raw[,3*i-1]=methy[,i]
  raw[,3*i]=dectP[,i]
  
}
colnames(raw)[seq(1,3*ncol,by=3)]<-paste(sam[match(colnames(beta),sam[,2]),1],"Unmethylated Signal",sep=" ")
colnames(raw)[seq(2,3*ncol,by=3)]<-paste(sam[match(colnames(beta),sam[,2]),1],"Methylated signal",sep=" ")
colnames(raw)[seq(3,3*ncol,by=3)]<-paste(sam[match(colnames(beta),sam[,2]),1],"Detection Pval",sep=" ")
rownames(raw)=rownames(beta)
write.table(raw,file="Matrix_signal_intensities.txt",sep="\t",quote=F,row.names=T,col.names=NA)






#
beta.v<-estimateBeta(idat.qtl, returnType=c("matrix"), offset = 100)
match(rownames(beta.v),)
BMIQ(beta.v,design.v,nL=3,doH=TRUE,nfit=5000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1)

sampleType <- pData(idat.qtl)$SampleType

allResult <- detectDMR.slideWin(idat.qtl, sampleType=sampleType)
allDMRInfo = identifySigDMR(allResult)
DMRInfo.ann <- annotateDMRInfo(allDMRInfo, "TxDb.Hsapiens.UCSC.hg19.knownGene")
export.DMRInfo(DMRInfo.ann, savePrefix="testExample")

library("xDb.Hsapiens.UCSC.hg19.knownGene")
phenoData<-c(rep(0,12),rep(1,12))
plotMethylationHeatmapByGene(selGene=6493,methyGenoSet=idat,phenoData=phenoData, includeGeneBody=TRUE,genomicFeature="TxDb.Hsapiens.UCSC.hg19.knownGene")



pdf("QC.before.pdf")
plotColorBias1D(idat)
plotColorBias2D(idat, selSample=1, cex=2)
boxplotColorBias(idat, channel='sum')
boxplot(estimateIntensity(idat))
dev.off()

pdf("QC.after.pdf")
plotColorBias1D(idat.cjc)
plotColorBias2D(idat.cjc, selSample=1, cex=2)
boxplotColorBias(idat.cjc, channel='sum')
boxplot(estimateIntensity(idat.cjc))
dev.off()


exprs=beta2m(myLoad$beta)
methylated=
unmethylated=
detection=myLoad$detP
methylated.N=
unmethylated.N=
assayData=
Champ2MethyLumiM<-new("MethyLumiM", exprs, methylated, unmethylated, detection, methylated.N, unmethylated.N,assayData)



showClass("MethyLumiM")



library("ChAMP")
library("outliers")
setwd("/home/sguo/dyh/idat")

setwd("G:\\dyh")
myLoad<-champ.load(directory = getwd(),filterBeads=TRUE,QCimages = F)
save(myLoad,file="myLoad.RData")
# turn off the plot if you run it on the server since of the problem of X11
myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = "sampleSheet.txt", resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = TRUE, filterXY = TRUE, QCimages = F, plotBMIQ = F)
save(myNorm,file="myNorm.RData")
TMR<-champ.TrueMethyl(beta.norm = myNorm$beta, pd = myLoad$pd, adjPVal = 0.5, adjust.method = "BH",compare.group = c("C", "T"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = TRUE)




load("G:\\geo\\GSE35069_matrix.Rdata")
data1 <- as.data.frame(exprs(GSE35069[[1]]))
phen <- pData(phenoData(GSE35069[[1]]))
phen1<-sapply(strsplit(as.character(phen$title),"_"),function(x) unlist(x)[1])

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\Annotation\\")
map<-read.table("GPL13534.sort.bed",sep="\t",as.is=T)

load("G:\\dyh\\myNorm.RData")
data2=myNorm$beta
phen2=c(rep("CA-RA-CD4+",12),rep("CA-HP-CD4+",12))

# merge GSE35069 and RA CD4+ array.

newdata<-data.frame(data1[match(rownames(data2),rownames(data1)),],data2)
data=t(na.omit(newdata))
pheno=c(phen1,phen2)
output="RA.CD4.pdf"
rlt<-PCAPlot(data,pheno,output,multifigure=T)


setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\RA\\analysis")
diff<-read.table("differential.analysis.txt",sep="\t",head=T,as.is=T)

library(gplots) 

mydata <- scale(t(data2))
sigcpg<-diff[which(diff$pvalue<0.005/nrow(diff)),1]
mydata1<-mydata[,na.omit(match(sigcpg,colnames(mydata)))]
which(rownames(mydata1)=="B4")
which(rownames(mydata1)=="B8")
mydata1<-mydata1[-c(4,8),]
rownames(mydata1)=c("RA0010","RA0002","RA0003","RA0005","RA0006","RA0007","RA0009","RA0001","RA0011","RA0012","HP0010","HP0002","HP0003","HP0004","HP0005","HP0006","HP0007","HP0008","HP0009","HP0001","HP0011","HP0012")

heatmap(mydata1,col=redgreen(10))

d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit,hang=-1,yaxt="n",main="",lwd=2.5,col="blue",cex=1.25) # display dendogram
axis(2,at=seq(0,100,by=20),lwd=2.5,cex=1.0)


ColorDendrogram(fit,y=c(rep(1,10),rep(2,12)),branchlength=11,labels=rownames(mydata))
ColorDendrogram

data=t(data2)
pheno=phen2


PCAPlot<-function(data,pheno,output,multifigure=T){
  pca <- prcomp(data,center=T,scale = F)  # Here, input file: row is individual and column is variable
  outputfile=paste(output,".pdf",sep="")
  pdf(outputfile)
  if(multifigure){
    par(mfrow=c(2,2),mar=c(4,4,4,4))  
  }
  plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
  axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
  var<-c()
  for(i in 1:length(pca$sdev)){var[i]<-sum((pca$sdev[1:i])^2)/sum((pca$sdev)^2)}
  plot(var,ylab="Total variance",xlab="Number of principle components",lwd=3,type="l",col="red")
  abline(h=0.8,col="grey",lty=2,lwd=2)
  abline(v=which(var>0.8)[1],col="grey",lty=2,lwd=2)
  scores <- data.frame(pheno, pca$x[,1:3])
  col = as.numeric(as.factor(pheno))
  plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=1.2,lwd=2)
  }
  legend("topleft",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",cex=1.05,pt.cex=1.2)
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n")
  dev.off()
}











