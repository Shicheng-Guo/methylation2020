RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}
library("Haplin")


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/GBM/GEO")
library("GEOquery")

GSE36194 <- getGEO("GSE36194")
data <- as.data.frame(exprs(GSE36194[[1]]))
phen <- pData(phenoData(GSE36194[[1]]))
GSE36194<-list()
GSE36194$beta=data
GSE36194$phen=phen
save(GSE36194,file="GSE36194.RData")

GSE53229 <- getGEO("GSE53229")
data <- as.data.frame(exprs(GSE53229[[1]]))
phen <- pData(phenoData(GSE53229[[1]]))
GSE53229<-list()
GSE53229$beta=data
GSE53229$phen=phen
save(GSE53229,file="GSE53229.RData")

GSE50923 <- getGEO("GSE50923")
data <- as.data.frame(exprs(GSE50923[[1]]))
phen <- pData(phenoData(GSE50923[[1]]))
GSE50923<-list()
GSE50923$beta=data
GSE50923$phen=phen
save(GSE50923,file="GSE50923.RData")


GSE46015 <- getGEO("GSE46015")
data <- as.data.frame(exprs(GSE46015[[1]]))
phen <- pData(phenoData(GSE46015[[1]]))
GSE46015<-list()
GSE46015$beta=data
GSE46015$phen=phen
save(GSE46015,file="GSE46015.RData")

GSE22867 <- getGEO("GSE22867")
data <- as.data.frame(exprs(GSE22867[[1]]))
phen <- pData(phenoData(GSE22867[[1]]))
GSE22867<-list()
GSE22867$beta=data
GSE22867$phen=phen
save(GSE22867,file="GSE22867.RData")


GSE119774 <- getGEO("GSE119774")
data <- as.data.frame(exprs(GSE119774[[1]]))
phen <- pData(phenoData(GSE119774[[1]]))
GSE119774<-list()
GSE119774$beta=data
GSE119774$phen=phen
save(GSE119774,file="GSE119774.RData")

GSE116298 <- getGEO("GSE116298")
data <- as.data.frame(exprs(GSE116298[[1]]))
phen <- pData(phenoData(GSE116298[[1]]))
GSE116298<-list()
GSE116298$beta=data
GSE116298$phen=phen
save(GSE116298,file="GSE116298.RData")


GSE36194 <- getGEO(filename="GSE36194_series_matrix.txt.gz",destdir="./")
data <- as.data.frame(exprs(GSE36194[[1]]))
phen <- pData(phenoData(GSE36194[[1]]))
GSE36194<-list()
GSE36194$beta=data
GSE36194$phen=phen
save(GSE36194,file="GSE36194.RData")

GSE53229<- getGEO(filename="GSE53229-GPL8490_series_matrix.txt.gz",destdir="./")
data <- as.data.frame(exprs(GSE116298[[1]]))
phen <- pData(phenoData(GSE116298[[1]]))
GSE116298<-list()
GSE116298$beta=data
GSE116298$phen=phen
save(GSE116298,file="GSE116298.RData")

GSE50923<- getGEO(filename="GSE50923_series_matrix.txt.gz",destdir="./")
data <- as.data.frame(exprs(GSE116298[[1]]))
phen <- pData(phenoData(GSE116298[[1]]))
GSE116298<-list()
GSE116298$beta=data
GSE116298$phen=phen
save(GSE116298,file="GSE116298.RData")


GSE46015<-read.table("GSE46015_non_normalized.txt",head=T,sep="\t")
GSE22867<-read.table("GSE22867_raw_illumina.txt",head=T,skip=4,sep="\t")



