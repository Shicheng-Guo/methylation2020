########################################################
########################################################
########################################################

library("impute")
library("DMRcate")
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

setwd("/home/sguo/tcga/lung")
load("TCGA.Methy450.LUAD.RData")
data<-RawNARemove(methdata$matrix)

rs<-grep("rs",rownames(data))
if(length(rs)>0){
data<-data[-rs,]  
}

myBetas<-impute.knn(data)$data
colnames(myBetas)<-colnames(data)
rownames(myBetas)<-rownames(data)

sampleid<-colnames(data)
idv<-unique(substr(sampleid,1,12))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(t1 %in% sampleid,t2 %in% sampleid)){
    pairidv<-c(pairidv,t1,t2)
  }
}
myBetas<-data[,match(pairidv,colnames(data))]
save(myBetas,file="TCGA.Methy450.LUAD.RM.Impute.RData")
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(substr(colnames(myMs.noSNPs),1,12))
type <- factor(substr(colnames(myMs.noSNPs),14,15))
design <- model.matrix(~patient + type)
myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",design=design, coef=ncol(design))
save(myannotation,file="TCGA.Methy450.LUAD.RM.Impute.myannotation.RData")

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
rlt<-dmrcoutput$results
write.table(rlt,file="TCGA.Methy450.LUAD.DMR.impute.txt",col.names=NA,row.names=T,sep="\t",quote=F)


load("TCGA.Methy450.LUAD.RM.Impute.RData")
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(substr(colnames(myMs.noSNPs),1,12))
type <- factor(substr(colnames(myMs.noSNPs),14,15))
load("TCGA.Methy450.LUAD.RM.Impute.myannotation.RData")

sigene<-c("CDO1","NUPR1","PTPN20B","HSPA2","IFFO1","ZNF597","NPTX2","ZNF132","PCDHB2","ZSCAN12",
  "C20orf151","CBLC","ZNF135","C3orf14","ZNF300","DNALI1","ZNF471","HIST1H2BH","CWH43",
  "HOXB4","CCDC106","PCDHB5")
ii<-na.omit(match(sigene,rlt[,1]))

for(i in 1:length(ii)){
filename=paste("luad",rlt[ii[i],1],"pdf",sep=".")
pdf(filename)
DMR.plot(dmrcoutput=dmrcoutput, dmr=ii[i], betas=myBetas,phen.col=as.numeric(type),pch=16, toscale=TRUE, plotmedians=TRUE,cex=0.8)
dev.off()
}

? DMR.plot

########################################################
########################################################
########################################################


library("impute")
library("DMRcate")
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

setwd("/home/sguo/tcga/lung")
load("TCGA.Methy450.LUSC.RData")
data<-RawNARemove(methdata$matrix)

rs<-grep("rs",rownames(data))
if(length(rs)>0){
  data<-data[-rs,]  
}

myBetas<-impute.knn(data)$data
colnames(myBetas)<-colnames(data)
rownames(myBetas)<-rownames(data)

sampleid<-colnames(data)
idv<-unique(substr(sampleid,1,12))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(t1 %in% sampleid,t2 %in% sampleid)){
    pairidv<-c(pairidv,t1,t2)
  }
}
myBetas<-data[,match(pairidv,colnames(data))]
save(myBetas,file="TCGA.Methy450.LUSC.RM.Impute.RData")
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(substr(colnames(myMs.noSNPs),1,12))
type <- factor(substr(colnames(myMs.noSNPs),14,15))
design <- model.matrix(~patient + type)
myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",design=design, coef=ncol(design))
save(myannotation,file="TCGA.Methy450.LUAD.RM.Impute.myannotation.RData")

load("TCGA.Methy450.LUAD.RM.Impute.myannotation.RData")
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
rlt<-dmrcoutput$results
write.table(rlt,file="TCGA.Methy450.LUSC.DMR.impute.txt",col.names=NA,row.names=T,sep="\t",quote=F)



load("TCGA.Methy450.LUSC.RM.Impute.RData")
load("TCGA.Methy450.LUSC.RM.Impute.myannotation.RData")

sigene<-c("CDO1","NUPR1","PTPN20B","HSPA2","IFFO1","ZNF597","NPTX2","ZNF132","PCDHB2","ZSCAN12",
          "C20orf151","CBLC","ZNF135","C3orf14","ZNF300","DNALI1","ZNF471","HIST1H2BH","CWH43",
          "HOXB4","CCDC106","PCDHB5")
ii<-na.omit(match(sigene,rlt[,1]))

match("IFFO1",rlt[,1])

for(i in 1:length(ii)){
  filename=paste("lusc",rlt[ii[i],1],"pdf",sep=".")
  pdf(filename)
  DMR.plot(dmrcoutput=dmrcoutput, dmr=ii[i], betas=myBetas,phen.col=as.numeric(type),pch=16, toscale=TRUE, plotmedians=TRUE,cex=0.8)
  dev.off()
}


