
setwd("/home/shg047/monod/oct/data")

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


Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

bed2cor<-function(bed){
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],"-",unlist(strsplit(x,"\t"))[3],sep="")})
  cor<-gsub("[ ]","",cor)
  return(cor)
}



file1<-read.table("WGBS_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",as.is=T,row.names=1,check.names=F)
file1<-data.matrix(file1[,13:56])  # remove H1 and WBC, cancer
colnames(file1)
colnames(file1)<-gsub("_","-",colnames(file1))
colname2<-unlist(lapply(colnames(file1),function(x) unlist(strsplit(x,"[.]"))[1]))
saminfo2<-read.table("/home/shg047/monod/phase2/newsaminfo.txt",head=T,sep="\t",as.is=T)
saminfo2<-saminfo2[na.omit(match(colname2,saminfo2[,1])),]
colnames(file1)<-saminfo2[,2]
colnames(file1)
tissue<-sort(c("Brain","Heart","muscle","Vessel","Kidney","Ovary","Esophagus","Lung","Liver","Pancreas","Stomach","Gastric","Intestine","Colon","Bladder"))
file1<-file1[,colnames(file1) %in% tissue]


# missing value distribution of file 1
pdf("WGBS.missing.value.distribution.pdf")
na_number<-unlist(apply(file1,1,function(x) sum(is.na(x))))
hist(na_number,xlim=c(0,50))
dev.off()



file2<-read.table("RRBS_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file2)
dim(file2)
# remove solid tissue
samplename1=sapply(strsplit(colnames(file2),"[.]"),function(x) unlist(x)[1])  # get sample name
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])        # get sample id
remove=c(samplename2[grep("6-T",samplename2)],samplename2[grep("PC-T",samplename2)],samplename2[grep("CTT-",samplename2)],samplename2[grep("7-T",samplename2)])
file2<-file2[,-match(remove,samplename2)]

samplename1=sapply(strsplit(colnames(file2),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
new<-read.table("/home/shg047/monod/phase2/saminfo.txt",sep="\t",as.is=T)
cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
samplename2<-gsub("6-P","CC-P",samplename2)
samplename2<-gsub("7-P","LC-P",samplename2)
samplename2<-gsub("6-T","CC-T",samplename2)
samplename2<-gsub("7-T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
colnames(file2)=samplename2

pdf("RRBS.missing.value.distribution.pdf")
na_number<-unlist(apply(file2,1,function(x) sum(is.na(x))))
hist(na_number)
dev.off()


# exactly same
cor1<-rownames(file1)
cor2<-rownames(file2)

file1<-file1[na.omit(match(cor2,cor1)),]
file2<-file2[cor2%in% cor1,]

data<-data.frame(file1,file2,check.names = F)

## optional
# colnames(data)[order(apply(data,2,function(x) sum(is.na(x))),decreasing=T)[1:25]]
# data<-data[,-order(apply(data,2,function(x) sum(is.na(x))),decreasing=T)[1:25]]

# f1<-c(grep("LC",colnames(data)))
# f1<-c(grep("CC",colnames(data)))
# f1<-c(grep("PC",colnames(data)))
f1<-c(grep("CC",colnames(data)),grep("LC",colnames(data)),grep("PC",colnames(data)))
f2<-c(grep("NC",colnames(data)))
f3<-(1:ncol(data))[-c(f1,f2)]
data=data[,c(f1,f3)]

data<-RawNARemove(data,missratio=0.4)
data<-impute.knn(data.matrix(data))$data

# optional
x1<-which(unlist(apply(data,1,function(x) sum(x>0.4)/ncol(data)))>0.7)
x2<-which(unlist(apply(data,1,function(x) sum(x<0.5)/ncol(data)))>0.7)
data<-data[-c(x1,x2),]

# esential
p.value<-unlist(apply(data,1,function(x) t.test(x[f1],x[f3])$p.value))
stastic<-unlist(apply(data,1,function(x) t.test(x[f1],x[f3])$statistic))

x1<-which(abs(stastic)>1.5)
x2<-which(p.value<0.005)
heatmapdata<-data[x1[x1 %in% x2],]
heatmapdata<-heatmapdata[,-c(36,43)]
# colnames(heatmapdata)[f1]<-"Cancer"
# colnames(heatmapdata)[f2]<-"Normal"

library("grDevices")
library("gplots")
filename=paste("Figure-cancer.signature-AK","pdf",sep=".")
pdf(filename)
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(heatmapdata,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.5,labRow=NA)
dev.off()



####
colnames(heatmapdata)
library("randomForest")

x=t(heatmapdata)
y=as.factor(colnames(heatmapdata))
rf<-randomForest(x=x,y=y,importance=T)

target<-rownames(rf$importance)[order(rf$importance[,4],decreasing = T)[1:20]]
x<-x[,match(target,colnames(x))]
rf<-randomForest(x=x,y=y,importance=T)
rf

library("grDevices")
library("gplots")
filename=paste("Figure-cancer.signature-2.rf","pdf",sep=".")
pdf(filename)
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(t(x),col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.5,labRow=NA)
dev.off()




heatmapdata2<-heatmapdata[order(rf$importance,decreasing=T)[1:20],]






setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\monod\\analysis\\phase3\\enrichmentAnalysis")
load("heatmap.data.RData")

data<-RawNARemove(data,missratio=0.4)
data<-impute.knn(data.matrix(data))$data

f1<-c(grep("LC",colnames(data)),grep("CC",colnames(data)),grep("PC",colnames(data)))
f2<-(1:ncol(data))[-f1]
p.value<-unlist(apply(data,1,function(x) t.test(x[f1],x[f2])$p.value))
stastic<-unlist(apply(data,1,function(x) t.test(x[f1],x[f2])$statistic))

x1<-which(abs(stastic)>1.5)
x2<-which(p.value<0.00005)


heatmapdata<-data[x1[x1 %in% x2],]
colnames(heatmapdata)[f1]<-"Cancer"
colnames(heatmapdata)[f2]<-"Normal"
dim(heatmapdata)

library("grDevices")
library("gplots")
filename=paste("Figure-cancer.signature","pdf",sep=".")
pdf(filename)
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(heatmapdata,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.2,labRow=NA)
dev.off()




