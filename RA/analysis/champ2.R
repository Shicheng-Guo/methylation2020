ManhattanmyDMP<-function(myDMP,output="ManhattanmyDMP.pdf"){
  library(qqman)
  SNP=rownames(myDMP)
  CHR=myDMP$CHR
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=myDMP$MAPINFO
  P=myDMP$P.Value
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=max(subset(myDMP,adj.P.Val<0.05)$P.Value)
  pdf(output)
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),lwd=2, suggestiveline=F)
  dev.off()
}

qqplot<-function(pvalues,output="qqplot.pdf"){
  library("Haplin")
  pdf(output)
  pQQ(pvalues, nlabs =length(pvalues), conf = 0.95) 
  dev.off()
}


library("ChAMP")
library("outliers")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")
champ.process(directory ="\\\\mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")
myLoad<-champ.load(directory = getwd(),filterBeads=TRUE,QCimages = F)
save(myLoad,file="myLoad.RData")
# turn off the plot if you run it on the server since of the problem of X11
myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = "sampleSheet.txt", resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = TRUE, filterXY = TRUE, QCimages = F, plotBMIQ = F)
save(myNorm,file="myNorm.RData")
TMR<-champ.TrueMethyl(beta.norm = myNorm$beta, pd = myLoad$pd, adjPVal = 0.5, adjust.method = "BH",compare.group = c("C", "T"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = TRUE)


library("GEOquery")
# GSE52826 <- getGEO("GSE52826",destdir="/home/sguo/monod/data/geo")
# save(GSE52826, file="GSE53045_matrix.Rdata")
load("GSE53045_matrix.Rdata")
data <- as.data.frame(exprs(GSE52826[[1]]))
phen <- pData(phenoData(GSE52826[[1]]))
# type1<-phen[match(colnames(data),rownames(phen)),]$characteristics_ch1
# PCAPlot(t(na.omit(data)),type1,output="esca.pca.GSE52826.pdf",legend.cex=0.5,multifigure=T)
# type2<-phen[match(colnames(data),rownames(phen)),]$characteristics_ch1.1

library("ChAMP")
library("doParallel")
setwd("~/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")
Dir="~/hpc/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat"

myLoad <- champ.load(Dir,filterBeads=TRUE,arraytype="EPIC")
pdf("MCaldwell.AMP.EPIC.QC.pdf")
champ.QC()
dev.off()

myNorm1 <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=6)
myNorm2 <- champ.refbase(beta=myNorm1,arraytype="EPIC")$CorrectedBeta

myDMP1 <- champ.DMP(beta = myNorm1,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 1)  # 2614
myDMP2 <- champ.DMP(beta = myNorm2,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 1)  # 2614

ManhattanmyDMP(myDMP1$Case_to_Control,output="manhattan.myDMP1.pdf")
qqplot(myDMP1$Case_to_Control$P.Value,output="qqplot.myDMP1.pdf") 
ManhattanmyDMP(myDMP2$Case_to_Control,output="manhattan.myDMP2.pdf")
qqplot(myDMP2$Case_to_Control$P.Value,output="qqplot.myDMP2.pdf") 

myDMR1 <- champ.DMR(beta=myNorm1,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=4,maxGap=1000)



write.table(myRefBase1$CellFraction,file="Table-S2.CellFraction.txt",sep="\t",quote=F,col.names = NA,row.names = T)
myNorm<-myRefBase$CorrectedBeta
myCellFraction<-champ.refbase(beta=myNorm,arraytype="EPIC")$CellFraction

pdf("MCaldwell.AMP.EPIC.SVD.pdf")
champ.SVD(beta=myNorm,pd=myLoad$pd)
dev.off()
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))

myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 0.05)  # 2614
ManhattanmyDMP(myDMP1$Case_to_Control,output="manhattan_raw.pdf")
qqplot(myDMP$Case_to_Control$P.Value,output="qqplot.pdf") 

myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=6)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$cmG1,arraytype="EPIC",adjPVal = 0.05)




write.table(myDMP$Old_Case_to_Old_Control,file="AtrialFibrillation.CaseControl.myDMP.OldCasevsOldControl.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myDMP$Young_Case_to_Young_Control,file="AtrialFibrillation.CaseControl.myDMP.YoungCasevsYoungControl.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myDMP$Old_Case_to_Young_Control,file="AtrialFibrillation.CaseControl.myDMP.Old_Case_to_Young_Control.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$cmG2,arraytype="EPIC",adjPVal = 0.1)
write.table(myDMP$Case_Before_to_Control_Before,file="AtrialFibrillation.CaseControl.myDMP.Case_Before_to_Control_Before.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMR1 <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=4,maxGap=1000)
write.table(myDMR$BumphunterDMR,file="AtrialFibrillation.CaseControl.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMR2 <- champ.DMR(beta=myNorm,pheno=myLoad$pd$cmG1,method="Bumphunter",compare.group=c("Old_Case","Old_Control"),arraytype="EPIC",minProbes=2,cores=4,maxGap=1000)
write.table(myDMR2$BumphunterDMR,file="AtrialFibrillation.myDMR.OldcaseVsOldControl.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMR3 <- champ.DMR(beta=myNorm,pheno=myLoad$pd$cmG1,method="Bumphunter",compare.group=c("Young_Case","Young_Control"),arraytype="EPIC",minProbes=2,cores=4,maxGap=1000)
write.table(myDMR3$BumphunterDMR,file="AtrialFibrillation.myDMR.Young_CaseVsYoung_Control.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMR3 <- champ.DMR(beta=myNorm,pheno=myLoad$pd$cmG2,method="Bumphunter",compare.group=c("Case_Before","Case_After"),arraytype="EPIC",minProbes=2,cores=5,maxGap=1000)
write.table(myDMR3$BumphunterDMR,file="AtrialFibrillation.myDMR.Case_Before_Case_After.txt",col.names = NA,row.names = T,quote=F,sep="\t")



write.table(myDMR,file="AtrialFibrillation.YoungOld.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
write.table(myDMR,file="AtrialFibrillation.Enrollment.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")



output<-methtools()


methtools<-function(){
output<-list()
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",adjPVal = 0.05)
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="EPIC",cores=6,maxGap=1000)

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="EPIC",minProbes=2,cores=6,maxGap=1000)
myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
ManhattanmyDMP(myDMP1$Case_to_Control,output="manhattan.pdf")
pQQ(myDMP$Case_to_Control$P.Value, nlabs =length(myDMP$Case_to_Control$P.Value), conf = 0.95, mark = F,output="qqplot.pdf")

output$myDMP<-myDMP
output$myDMR<-myDMR
output$myEpiMod<-myEpiMod
output$myBlock<-myBlock
output$myGSEA<-myGSEA
output$myCNA<-myCNA
return(output)
}







write.table(myBlock$Block,file="AtrialFibrillation.CaseControl.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.YoungOld.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.Enrollment.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
write.table(myebayGSEA,file="AtrialFibrillation.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase1 <- champ.refbase(beta=myNorm,arraytype="EPIC")
myRefBase2 <- champ.refbase(beta=myNorm,arraytype="450K")


BiocManager::install("openxlsx") 
BiocManager::install("xlsx") 
BiocManager::install("ggpubr") 

library("openxlsx")
library("xlsx")

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Ingrid/MCaldwell-Sept27-17-HuMethEPIC/Raw_Data/idat")
data<-openxlsx::read.xlsx("Tables.xlsx",sheet=1)

pdf("MethylationAge.pdf")
boxplot(DNAmAge_Norm~cmG1,data=data,col=c(2,3,2,3),xlab="Group",ylab="DNA Methylation Age")
boxplot(DNAmAge_Norm~cmG2,data=subset(data,Sample_Group=="Case"),col=c(2,3),xlab="Group",ylab="DNA Methylation Age")
boxplot(DNAmAge_Norm~Sample_Group,data=data,col=2:5,xlab="Group",ylab="DNA Methylation Age")
dev.off()

saminfo<-openxlsx::read.xlsx("Tables.xlsx",sheet=1)
data<-openxlsx::read.xlsx("Tables.xlsx",sheet=3,rowNames=T)
rownames(data)<-paste(saminfo[match(rownames(data),saminfo$Sample_Name),]$Sample_Group,rownames(data),sep="-")
colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
heatmap.2(data.matrix(data),trace="none",density.info="none",keysize=0.9,col=colors)

cc<-saminfo[match(rownames(data),saminfo$Sample_Name),]$Sample_Group
P<-apply(data,2,function(x) t.test(x~cc))



kruskal.test(DNAmAge_Norm ~ cmG1, data = data)
kruskal.test(DNAmAge_Norm ~ cmG2, data = data)
aov(DNAmAge_Norm ~ cmG1, data = data)
ggpubr::ggboxplot(my_data, x = "group", y = "weight", color = "group")








