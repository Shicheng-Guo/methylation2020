#####################################################################################
###   Title : Methylation Block Region analysis based on Methylation 450K beadchip
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Section 1.  Collect Methylation 450K array and Different methylation loci
##   Section 2.   Collect Gene expression and Different expression genes
###   Section 3.  Correlation between gene expression and methylation
###   Section 4.  Summary and Plot 
#####################################################################################

library("impute")
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

map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")

setwd("/home/sguo/tcga/esca/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3")   # need update
file=list.files(pattern="jhu-usc.edu_ESCA.HumanMethy*")                             # need update
sampleinfo=read.table("/home/sguo/tcga/esca/file_manifest.txt",head=T,sep="\t")   # need update
sampleid=sampleinfo[match(file,sampleinfo[,ncol(sampleinfo)]),5]
sampleid

data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data<-cbind(data,tmp[,2])
  print(i)
}
rownames(data)=tmp[,1]
colnames(data)=sampleid

idv<-unique(substr(sampleid,1,12))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(t1 %in% sampleid,t2 %in% sampleid)){
    pairidv<-c(pairidv,t1,t2)
  }
}
methdata<-data[,match(pairidv,colnames(data))]
save(methdata,file="/home/sguo/tcga/esca/TCGA.Methy450.ESCA.RData")
methdata=RawNARemove(methdata)
newmethdata<-impute.knn(methdata)$data
rownames(newmethdata)=rownames(methdata)
colnames(newmethdata)=colnames(methdata)
save(newmethdata,file="/home/sguo/tcga/esca/TCGA.Methy450.ESCA.RM.Impute.RData")
type<-substr(colnames(newmethdata),14,15)
table(type)
data=newmethdata
data<-data+matrix(rnorm(length(data),0.0001,0.0001),dim(data)[1],dim(data)[2])   # row=gene, col=inv
output<-matrix(NA,dim(data)[1],5)   # set output matrix ()
x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive
x2<-which(type==names(table(type))[2])   # type 2, normal or resistant
for(i in 1:dim(data)[1]){
  if(! any(all(is.na(data[i,x1])),all(is.na(data[i,x2])))){ 
    tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=T, na.action=na.omit))
    output[i,1]<-tmp1$p.value
    output[i,2]<-mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2]))
    output[i,3]<-"ttest"
    output[i,4]<-mean(as.numeric(data[i,x1]))
    output[i,5]<-mean(as.numeric(data[i,x2]))
    print(i)
  }
}
out<-cbind(rownames(data),output)
colnames(out)<-c("GeneName","Pvalue","Statistic","test","mean1","mean2")
write.table(out,file="/home/sguo/tcga/esca/TCGA.Methy450.diff.esca.txt",sep="\t",quote=F,col.names=NA,row.names=T)
p<-as.numeric(out[,2])
rlt.sig<-out[which(p<0.05/length(p)),]

write.table(rlt.sig,file="/home/sguo/tcga/esca/TCGA.Methy450.sig.diff.esca.txt",sep="\t",quote=F,col.names=NA,row.names=T)  # update

#################################################################################
########################  Gene expression section  ##############################
#################################################################################

setwd("/home/sguo/tcga/esca/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3")
file=list.files(pattern="*rsem.genes.normalized_results")
sampleinfo=read.table("/home/sguo/tcga/esca/file_manifest.txt",head=T,sep="\t")
sampleid=sampleinfo[match(file,sampleinfo[,ncol(sampleinfo)]),5]

data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=T,sep="\t",as.is=F)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data<-cbind(data,tmp[,2])
  print(i)
}
rownames(data)=tmp[,1]
colnames(data)=sampleid
idv<-unique(substr(sampleid,1,12))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(t1 %in% sampleid,t2 %in% sampleid)){
    pairidv<-c(pairidv,t1,t2)
  }
}
expdata<-data[,match(pairidv,colnames(data))]
save(expdata,file="/home/sguo/tcga/esca/TCGA.RNAseqV2.ESCA.RData")

expdata=RawNARemove(expdata)
newexpdata<-impute.knn(expdata)$data
rownames(newexpdata)=rownames(expdata)
colnames(newexpdata)=colnames(expdata)
save(newexpdata,file="/home/sguo/tcga/esca/TCGA.RNAseqV2.ESCA.RM.Impute.RData")
type<-substr(colnames(newexpdata),14,15)
table(type)
data=newexpdata
# add noise 
data<-data+abs(matrix(rnorm(length(data),0.001,0.001),dim(data)[1],dim(data)[2]))   # row=gene, col=inv
# obtain x,y index for non-pair t or wilcox test
data<-log(data,2)
output<-matrix(NA,dim(data)[1],5)   # set output matrix ()
x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive
x2<-which(type==names(table(type))[2])   # type 2, normal or resistant
for(i in 1:dim(data)[1]){
  if(! any(all(is.na(data[i,x1])),all(is.na(data[i,x2])))){ 
    tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=T, na.action=na.omit))
    output[i,1]<-tmp1$p.value
    output[i,2]<-mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2]))
    output[i,3]<-"ttest"
    output[i,4]<-mean(as.numeric(data[i,x1]))
    output[i,5]<-mean(as.numeric(data[i,x2]))
    print(i)
  }
}
out<-cbind(rownames(data),output)
colnames(out)<-c("GeneName","Pvalue","Statistic","test","mean1","mean2")
write.table(out,file="/home/sguo/tcga/esca/TCGA.RNAseqV2.diff.esca.txt",sep="\t",quote=F,col.names=NA,row.names=T)  # update
p=as.numeric(out[,2])
rlt.sig<-out[which(p<0.05/length(p)),]
write.table(rlt.sig,file="/home/sguo/tcga/esca/TCGA.RNAseqV2.sig.diff.esca.txt",sep="\t",quote=F,col.names=NA,row.names=T)  # update

################################################################################################
########################  Methylation Gene expression Regulation  ##############################
################################################################################################
map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")
load("TCGA.Methy450.ESCA.RM.Impute.RData")
load("TCGA.RNAseqV2.ESCA.RM.Impute.RData")
methydata<-newmethdata
expdata<-newexpdata[,na.omit(match(colnames(methydata),colnames(newexpdata)))]
methydata<-methydata[,colnames(methydata) %in% colnames(expdata)]

expdata[1:5,1:5]
methydata[1:5,1:5]
dim(expdata)
dim(methydata)

gene1<-map[match(rownames(methydata),map[,4]),5]
gene2<-sapply(rownames(expdata),function(x) unlist(strsplit(x,"[|]"))[1])
gene<-unique(as.character(gene2[which(gene2 %in% gene1)]))

genename<-c()
cpgname<-c()
beta<-c()
pvalue<-c()

for (i in 1:length(gene)){
  nro<-which(gene1==gene[i])
  for (j in 1:length(nro)){
    pvaluetmp<-try(summary(lm(expdata[match(gene[i],gene2),]~methydata[nro[j],]))$coefficients[2,4],silent = TRUE)
    betatmp<-try(summary(lm(expdata[match(gene[i],gene2),]~methydata[nro[j],]))$coefficients[2,1],silent = TRUE)
    if(! is.na(as.numeric(pvaluetmp))){
      genename<-c(genename,gene[i])
      cpgname<-c(cpgname,rownames(methydata)[nro[j]])
      beta<-c(beta,betatmp)
      pvalue<-c(pvalue,pvaluetmp)
    }
  }
  print(i)
}
rlt<-data.frame(genename,cpgname,beta,pvalue)
write.table(rlt,file="TCGA.ESCA.methylation-expression-regulation.esca.txt",sep="\t",quote=F,col.names=T,row.names=F)
rlt.sig<-rlt[rlt$pvalue<0.05/(nrow(rlt)),]
write.table(rlt.sig,file="TCGA.ESCA.Sig.methylation-expression-regulation.txt",sep="\t",quote=F,col.names=T,row.names=F)


################################################################################################
########################  Methylation Gene expression Regulation  ##############################
################################################################################################
map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")

sig1<-read.table("TCGA.ESCA.methylation-expression-regulation.luad.txt",sep="\t",head=T)
sig2<-read.table("tcga.Methy450.sig.diff.esca.txt",sep="\t",head=T,row.names=1)
sig3<-read.table("tcga.RNAseqV2.sig.diff.esca.txt",sep="\t",head=T,row.names=1)

sig2<-data.frame(sig2,gene=map[na.omit(match(sig2[,1],map[,4])),1:7])
write.table(sig2,file="/home/sguo/tcga/esca/TCGA.Methy450.diff.esca.txt",sep="\t",quote=F,col.names=NA,row.names=T)

head(sig1)
head(sig2)
head(sig3)

match(sig2[,1],sapply(as.character(sig3[,1]),function(x) unlist(strsplit(x,"[|]"))[1]))


