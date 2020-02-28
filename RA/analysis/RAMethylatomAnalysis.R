##########################################################################################
###   Title : Genome-wide DNA methylation analysis for RA patient with ChAMP package
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Section 1. Function predifinition 
###   Section 2. Data Cleaning
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
##########################################################################################
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\RA\\analysis")
sam<-read.table("sam.txt",head=T,sep="\t",as.is=T)
head(sam)
setwd("G:/dyh/")
library("ChAMP")
library("outliers")
setwd("/home/sguo/dyh/idat")
myLoad<-champ.load(directory = getwd(),filterBeads=TRUE,QCimages = F)
save(myLoad,file="myLoad.RData")
# turn off the plot if you run it on the server since of the problem of X11
case1<-myLoad$beta[,1:12]
con1<-myLoad$beta[,13:24]
myLoad$beta[,1:12]<-apply(case1,1,function(x) rm.outlier(x, fill = T, median = T, opposite = FALSE))
myLoad$beta[,13:24]<-apply(con1,1,function(x) rm.outlier(x, fill = T, median = T, opposite = FALSE))
myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = "sampleSheet.txt", resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = TRUE, filterXY = TRUE, QCimages = F, plotBMIQ = F)
save(myNorm,file="myNorm.RData")
TMR<-champ.TrueMethyl(beta.norm = myNorm$beta, pd = myLoad$pd, adjPVal = 0.05, adjust.method = "BH",compare.group = c("C", "T"), resultsDir = paste(getwd(), "resultsChamp2", sep = "/"),bedFile = TRUE)

library("ChAMP")
library("outliers")
setwd("/home/sguo/dyh/idat")
load("myLoad.RData")
load("myNorm.RData")
TMR<-champ.TrueMethyl(beta.norm = Normbeta, pd = myLoad$pd, adjPVal = 0.05, adjust.method = "BH",compare.group = c("C", "T"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = TRUE)


champ.CNA(batchCorrect=FALSE,sampleCNA=T,groupFreqPlots=T)
champ.process(directory = getwd(),adjPVal = 0.5,QCimages = F,,plotSample =F,groupFreqPlots=F)

x<-c(0.1,0.12,0.15,0.3,0.2,0.5,0.9,0.14,0.3,0.4,0.3,0.4)
t<-chisq.out.test(x, variance=var(x), opposite = FALSE)



source("http://bioconductor.org/biocLite.R")
library("heatmap.plus")
library("gplots")

setwd("/home/sguo/mice")
raw.sig<-read.table("matrix.detail.mice.alicey.txt",head=T,sep="\t")
data<-data.matrix(raw.sig[,c(7,4,8,9)])
clab<-c("x","129_mESC","B6_mESC","miPS_1","miPS_2","miPS_2","miPS_B3","SCNT_B12","SCNT_NB3","SCNT_P7C","SCNT_P8B")
colnames(data)=c("miPS_1","miPS_2","SCNT_1","SCNT_2")
pdf("mice.sig.pdf")
par(mar=c(1,1,1,1))
heatmap.plus(data,cex.main=0.2,cexRow=0.2,cexCol=0.75,col=redgreen(20));
dev.off()





