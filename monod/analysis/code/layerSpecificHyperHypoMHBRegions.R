########################################################################################
###   Title: layer specific methylation haplotype
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 11/9/2015
########################################################################################

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

###################################################################################################################

setwd("/home/shg047/monod/dec")
infile="WGBS_methHap_load_matrix_20Oct2015.txt";
file1<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)

# miss value detection and imputation
library("impute")
f2<-RawNARemove(file1,missratio=0.3)
f2<-impute.knn(data.matrix(f2))$data

colnames(f2)
library("preprocessCore")
f2.t1<-normalize.quantiles(f2[,13:58])
library("sva")
batch=c(rep(1,10),rep(2,36))
f2.t2<-ComBat(f2.t1, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)
f2[,13:58]<-f2.t2


# re-assign colnames
colnames(f2) 
colnames(f2)<-gsub("_","-",colnames(f2))
colname2<-unlist(lapply(colnames(f2),function(x) unlist(strsplit(x,"[.]"))[1]))
colname2
colnames(f2)<-colname2
# be sure all the sample information has been stored in the following database
saminfo2<-read.table("/home/shg047/monod/phase2/newsaminfo.txt",head=T,sep="\t",as.is=T)
saminfo2<-saminfo2[match(colname2,saminfo2[,1]),]
saminfo2
colnames(f2)<-saminfo2[,2]

saminfo3<-read.table("/home/shg047/monod/saminfo/tissue2Layer.txt",head=T,sep="\t",as.is=T)
f2<-f2[,saminfo2[,2] %in% saminfo3[,1]]
fn<-f2
colnames(fn)<-saminfo3[match(colnames(fn),saminfo3[,1]),2]
group=names(table(colnames(fn)))
index=colnames(fn)
gsi<-c()
gmaxgroup<-c()
pvalue=apply(fn,1,function(x) summary(aov(x~index))[[1]][["Pr(>F)"]][1])

pvalue=apply(fn,1,function(x) summary(aov(x~index))[[1]][[5]][1])

SigDiffMHBANOVA<-fn[match(names(which(pvalue<9.223561e-07)),rownames(fn)),]
save(SigDiffMHBANOVA,file="SigDiffMHBANOVA.RData")

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\layer_specfic_mhl\\anova")
library("gplots")
load("SigDiffMHBANOVA.RData")
SigDiffMHBANOVA[SigDiffMHBANOVA<0]<-0
SigDiffMHBANOVA[SigDiffMHBANOVA>1]<-1

SigDiffMHBANOVA<-SigDiffMHBANOVA[,order(colnames(SigDiffMHBANOVA))]
pdf("Figure.supervised.layer.mhl.single.cpg.heatmap.analysis.combat.quantile.pdf")
col=colorRampPalette(c("yellow", "blue"))(20) 
rlt<-heatmap.2(data.matrix(SigDiffMHBANOVA),col=col,distfun=distFun,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.85,labRow=NA)
dev.off()

? heatmap.2

rlt=data.frame(region=rownames(fn),group=gmaxgroup,GSI=gsi)
write.table(rlt,file="Table.GSI.layer.mhl.WGBS.Remove.H1.WBC.rlt.txt",col.names=T,row.names=F,quote=F,sep="\t")






