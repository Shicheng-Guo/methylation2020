#!/usr/bin/env Rscript
# please copy monod package from: /media/Home_Raid1/shg047/monod_1.1.tar.gz and install.packages("monod_1.1.tar.gz")

library("monod")
library("ggplot2")
library("reshape2")
#! Remember all the advance analysis should be conducted in Genome-miner so that everything can be documented.
# system("scp /oasis/tscc/scratch/shg047/monod/hapinfo/N37Salk.saminfo shg047@genome-miner.ucsd.edu:/media/NAS1/shg047/monod/hapinfo")
# system("scp /media/NAS1/shg047/tcga/pancancer/colon.450h.RData shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/monod/hapinfo")
#args = commandArgs(trailingOnly=TRUE)

gsi<-function (matrix){
  group = names(table(colnames(matrix)))
  index = colnames(matrix)
  GSI <- c()
  gmaxgroup <- c()
  refMean<-c()
  refMax<-c()
  for (i in 1:nrow(matrix)) {
    gsit <- 0    
    refmean<-tapply(as.numeric(matrix[i, ]), index, function(x) mean(x, na.rm = T))
    refmax<-refmean[which.max(refmean)]
    gmax <- names(refmax)
    for (j in 1:length(group)) {
      tmp <- (1 - 10^(mean(na.omit(as.numeric(matrix[i, which(index == group[j])])), na.rm = T))/10^(mean(na.omit(as.numeric(matrix[i,which(index == gmax)])))))/(length(group) - 1)
      gsit <- c(gsit, tmp)
    }
    if(sum(refmean>0.3,na.rm=T)>1){
      gmax<-gsub(" ","",paste(unique(c(gmax,names(refmean)[which(refmean>0.3)])),',',collapse =""))
    }
    gmaxgroup <- c(gmaxgroup, gmax)
    GSI <- c(GSI, sum(gsit, na.rm = T))
    refMean<-c(refMean,gsub(" ","",paste(round(refmean,5),',',collapse ="")))
    refMax<-c(refMax,refmax)    
  }  
  rlt = data.frame(region = rownames(matrix), group = gmaxgroup, GSI = GSI, refMax=refMax,refMean=refMean)
  return(rlt)
}

rename<-function (data){
  Data = data[, grep("STL|N37|ENC|SRX|age|new|centenarian|WB|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P", 
                     colnames(data))]
  colnames(Data)[grep(".", colnames(Data))] <- unlist(lapply(colnames(Data)[grep(".", 
                                                                                 colnames(Data))], function(x) unlist(strsplit(x, ".hapInfo|.sorted"))[1]))
  colnames(Data) <- gsub("[.]", "-", colnames(Data))
  colnames(Data)[grep("age|new|centenarian|middle", colnames(Data))] <- "WBC"
  colnames(Data)[grep("X6.T|CTT|CCT", colnames(Data))] <- "CCT"
  colnames(Data)[grep("X7.T|LCT", colnames(Data))] <- "LCT"
  colnames(Data)[grep("N37|STL|ENC", colnames(Data))] <- as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC", 
                                                                                                        colnames(Data))], saminfo[, 1]), 2])
  colnames(Data)[grep("X6.P|RRBS.6P", colnames(Data))] <- "CCP"
  colnames(Data)[grep("X7.P|RRBS.7P", colnames(Data))] <- "LCP"
  colnames(Data)[grep("NC.P", colnames(Data))] <- "NP"
  colnames(Data)[grep("X7.T|SRX381722|SRX381719|SRX381716", 
                      colnames(Data))] <- "LCT"
  colnames(Data)[grep("X6.T|CTT|SRX381569", colnames(Data))] <- "CCT"
  return(Data)
}

colon2hyphen<-function(x){
  colon<-lapply(x,function(x) unlist(strsplit(x,":")))
  cor<-lapply(colon,function(x) paste(x[1],":",x[2],"-",as.numeric(x[3])-1,sep=""))
  return(unlist(cor))
}

gsi2bio<-function(gsirlt){
  cor2bed<-function(cor){
    cor<-as.character(cor)
    a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
    bed<-matrix(a,ncol=3,byrow=T)
    return(data.frame(bed))
  }
  bio<-data.frame(cor2bed(gsirlt[,1]),gsirlt[,2:5])
  rownames(bio)<-gsirlt[,1]
  return(bio)
}

Zscore<-function(ccpr,npr){
  Zmax<-matrix(0,nrow=nrow(ccpr),ncol=ncol(ccpr))
  for(i in 1:ncol(ccpr)){
    for(k in 1:nrow(npr)){
      idx<-1:ncol(npr)
      zmp <- (ccpr[k,i] - mean(npr[k, idx]))/(sd(npr[k, idx])*sqrt((length(npr[k, idx])-1)/(length(npr[k, idx]))))
      Zmax[k,i]=zmp 
    }
  }
  rownames(Zmax)=rownames(ccpr)
  colnames(Zmax)=colnames(ccpr)
  Zmax
}


ZscorePredictionPlasma<-function(test,ncp,bio,tt=0.01){
  rlt<-list()
  test<-test[rownames(test) %in% rownames(bio),]
  ncp<-ncp[rownames(ncp) %in% rownames(bio),] 
  ccp<-test[rownames(test) %in% rownames(ncp),]
  npr<-scale(apply(ncp,2,function(x) table(unlist(lapply(bio[match(rownames(ncp)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,",")))))))
  ccpr<-scale(apply(ccp,2,function(x) table(unlist(lapply(bio[match(rownames(ccp)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,",")))))))
  rlt$nprscore<-Zscore(npr,npr)
  rlt$ccprscore<-Zscore(ccpr,npr)
  rlt$ccpvalue<-pnorm(-rlt$ccprscore)
  rlt$npvalue<-pnorm(-rlt$nprscore)
  rownames(rlt$ccpvalue)<-rownames(rlt$ccpvalue)
  colnames(rlt$ccpvalue)<-colnames(rlt$ccpvalue)
  rownames(rlt$npvalue)<-rownames(rlt$npvalue)
  colnames(rlt$npvalue)<-colnames(rlt$npvalue)
  rlt1<-apply(Zscore(ccpr,npr),2,function(x) rownames(npr)[which.max(x)])
  rlt2<-apply(Zscore(npr,npr),2,function(x) rownames(npr)[which.max(x)])
  rlt$test<-rlt1
  rlt$background<-rlt2
  return(rlt)
}

ZscorebaseTest<-function(data,bio){
  Sum<-c()
  for(tt in seq(0,0.6,by=0.01)){
    input<-data[match(rownames(bio),rownames(data)),]
    set.seed(1)
    np<-input[,grep("NC.P",colnames(input))]  
    ccp<-input[,grep(".6P|X6.P",colnames(input))]  
    lcp<-input[,grep(".7P|X7.P",colnames(input))] 
    train<-sample(1:75,45)
    trainset<-np[,train]
    testset<-np[,(1:75)[-train]]
    npr<-apply(trainset,2,function(x) table(unlist(lapply(bio[match(rownames(trainset)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,","))))))
    ccpr<-apply(ccp,2,function(x) table(unlist(lapply(bio[match(rownames(ccp)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,","))))))
    lcpr<-apply(lcp,2,function(x) table(unlist(lapply(bio[match(rownames(ccp)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,","))))))
    nnpr<-apply(trainset,2,function(x) table(unlist(lapply(bio[match(rownames(trainset)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,","))))))
    rlt1<-sum(apply(Zscore(ccpr,npr),2,function(x) which.max(x)==2))/ncol(ccp)
    rlt2<-sum(apply(Zscore(lcpr,npr),2,function(x) which.max(x)==6))/ncol(lcp)
    rlt3<-sum(apply(Zscore(nnpr,npr),2,function(x) which.max(x)==10))/ncol(testset)
    Sum<-rbind(Sum,c(ccp=rlt1,lcp=rlt2,npc=rlt3,total=rlt1+rlt2+rlt3))
  }
  return(Sum)
}

setwd("/media/NAS1/shg047/monod/hapinfo")
#system("scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/saminfo/N37Salk.saminfo ./")
#system("scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/hapinfo/MHL4.RData ./")
#data1<-read.table("/media/Home_Raid1/shg047/NAS1/monod/mhl/mhl.mhbs",head=T,row.names=1)
#rownames(data1)<-colon2hyphen(rownames(data1))
load("./MHL4.RData")
colnames(data)<-gsub(".sorted.clipped.bam.hapInfo.txt.hap","",colnames(data))
colnames(data)<-gsub(".hapInfo.txt.hap","",colnames(data))
saminfo<-read.table("/media/NAS1/shg047/monod/hapinfo/N37Salk.saminfo",sep="\t")
load("gsirlt.RData")              # selected biomarker for each reference
nrow(gsirlt)
gsirlt1<-subset(gsirlt,refMax>0.6 & GSI>0.6)
gsirlt1$group<-as.character(gsirlt1$group)
nrow(gsirlt1)
bio<-gsi2bio(gsirlt1)

np<-data[,grep("NC.P",colnames(data))]  
train<-sample(1:75,45)
ncp<-np[,train]
nncp<-np[,(1:75)[-train]]
ccp<-data[,grep(".6P|X6.P",colnames(data))]  
lcp<-data[,grep(".7P|X7.P",colnames(data))] 

rlt1<-ZscorePredictionPlasma(ccp,ncp,bio,tt=0.01)
rlt2<-ZscorePredictionPlasma(lcp,ncp,bio,tt=0.01)
rlt3<-ZscorePredictionPlasma(nncp,ncp,bio,tt=0.01)

acc1<-unlist(lapply(apply(rlt1$ccpvalue,2,function(x) rownames(rlt1$ccpvalue)[which(x<0.005)]),function(x) grep("Colon|Intestine",x)[1]))
acc2<-unlist(lapply(apply(rlt2$ccpvalue,2,function(x) rownames(rlt1$ccpvalue)[which(x<0.005)]),function(x) grep("Lung",x)[1]))
acc3<-unlist(lapply(apply(rlt3$ccpvalue,2,function(x) rownames(rlt1$ccpvalue)[which(x<0.005)]),function(x) grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach",x)[1]))
acc1
acc2
acc3

