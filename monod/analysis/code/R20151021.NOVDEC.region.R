########################################################################################
###   Title: Group Specificity Index (GSI) for Genome-wide Methylation Haplotype dataset
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 9/1/2015
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
bedwithgap<-function(bed,gap){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
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
  plot(var,ylab="total variance",xlab="number of principle components",lwd=2,type="l")
  abline(h=0.8,col="grey",lty=2)
  abline(v=which(var>0.8)[1],col="grey",lty=2)
  scores <- data.frame(pheno, pca$x[,1:3])
  col = as.numeric(as.factor(pheno))
  plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.8,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n")
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n")
  dev.off()
}

#########################################################################################
setwd("/media/Home_Raid1/dinh/WGBS_HaploInfo/scripts/Sharing/")
infile1="WGBS_methHap_load_matrix_16Oct2015.txt";
file1<-read.table(infile1,head=T,sep="\t",row.names=1,as.is=T,check.names=F)
infile2="WGBS_methHap_load_matrix_20Oct2015.txt";
file2<-read.table(infile2,head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file1<-file1[na.omit(match(rownames(file2),rownames(file1))),]
file2<-file2[rownames(file2) %in% rownames(file1),]
sum(is.na(file1))/(nrow(file1)*ncol(file1))
sum(is.na(file2))/(nrow(file1)*ncol(file2))
colnames(file1) %in% colnames(file2)
file2<-file2[,match(colnames(file1),colnames(file2))]
colnames(file1) == colnames(file2)
# compare two mhl matrix based on length>=1 and >=3
cor<-c()
for(i in 1:ncol(file1)){
  tmp<-cor(file1[,i],file2[,i],use="complete.obs")
  cor<-c(cor,tmp)
}
test<-data.frame(file1[,1],file2[,1])
colnames(test)<-c("Colon_L>=3","Colon_L>=1")
rownames(test)<-rownames(file1)
test[1:40,]
test<-data.frame(file1[,60],file2[,60])
colnames(test)<-c("Colon_L>=3","Colon_L>=1")
rownames(test)<-rownames(file1)
test[1:40,]
pdf("smoothScatter.pdf")
smoothScatter(test)
dev.off()
pdf("correlation.betweenL3L1.pdf")
barplot(cor,ylim=c(0,1.1),col="blue")
dev.off()

#########################################################################################
# average methylation level analysis

setwd("/home/shg047/monod/dec")
infile="WGBS_aveMeth_load_matrix_20Oct2015.txt";
file1<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)

# miss value detection and imputation
library("impute")
f2<-RawNARemove(file1,missratio=0.05)
f2<-impute.knn(data.matrix(f2))$data

sum(is.na(file1))/(nrow(file1)*ncol(file1))
sum(is.na(f2))/(nrow(f2)*ncol(f2))

library("preprocessCore")
f2[,13:61]<-normalize.quantiles(f2[,13:61])

f3 <- t(na.omit(f2)) # listwise deletion of missing
# f3 <- scale(f3) # standardize variables
# colnames(f3)=colnames(file1)
# colnames(f3)=colnames(f2)

# check the dendrogram after missing, quantile
histFfile=paste("Figure.Cluster.aml.Dendrogram.ward.D.plot",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
pdf(histFfile)
plot(fit,cex=0.5) # display dendogram
dev.off()

histFfile=paste("Figure.Cluster.Dendrogram.aml.plot.average",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()


#########################################################################################
setwd("/home/shg047/monod/dec")
infile="WGBS_methHap_load_matrix_20Oct2015.txt";
file1<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)

# miss value detection and imputation
library("impute")
f2<-RawNARemove(file1,missratio=0.3)
f2<-impute.knn(data.matrix(f2))$data

f2.t1<-normalize.quantiles(f2[,13:58])
library("sva")
batch=c(rep(1,10),rep(2,36))
f2.t2<-ComBat(f2.t1, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)
f2[,13:58]<-f2.t2

# re-plot the cluster by sequencing sample Tag
colnames(f2)<-colnames(file1)
f3 <- t(na.omit(f2)) # listwise deletion of missing
histFfile=paste("Figure.Cluster.mhl.Dendrogram.combat.seqTag.ward.D.plot",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()
histFfile=paste("Figure.Cluster.mhl.Dendrogram.combat.seqTag.plot.average",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()

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
colnames(f2)

# check the dendrogram after missing, quantile
f3 <- t(na.omit(f2)) # listwise deletion of missing
# f3 <- scale(f3) # standardize variables
sum(is.na(file1))/(nrow(file1)*ncol(file1))
sum(is.na(f2))/(nrow(f2)*ncol(f2))
histFfile=paste("Figure.hist.combat",infile,"pdf",sep=".")
pdf(histFfile)
boxplot(f2,outline=F,horizontal=T,notch=F,las=1,cex.axis=0.65,col="blue")
dev.off()
histFfile=paste("Figure.Cluster.mhl.Dendrogram.combat.ward.D.plot",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()
histFfile=paste("Figure.Cluster.mhl.Dendrogram.combat.plot.average",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()


# tissue specific hml
f2<-f2[,13:58]
group=names(table(colnames(f2)))
index=colnames(f2)
gsi<-c()
gmaxgroup<-c()
for(i in 1:nrow(f2)){
  gsit<-0
  gmax<-names(which.max(tapply(as.numeric(f2[i,]),index,mean)))
  for(j in 1:length(group)){
    tmp<-(1-10^(mean(f2[i,][which(index==group[j])]))/10^(mean(f2[i,][which(index==gmax)])))/(length(group)-1)
    gsit<-gsit+tmp
  }
  gmaxgroup<-c(gmaxgroup,gmax)
  gsi<-c(gsi,gsit)
  print(c(gmax,gsit))
}
rlt=data.frame(region=rownames(f2),group=gmaxgroup,GSI=gsi)

sum(gsi>0.6)
write.table(rlt,file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",col.names=T,row.names=F,quote=F,sep="\t")


# Figure Distribution of GSI (without H1,WBC,Cancer Cells)
data<-read.table(file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
pdf("Figure.GSI.hist.distribution.pdf")
hist(data$GSI,breaks=100,xlim=c(0,0.9),ylim=c(0,2000),col="red")
dev.off()

# Table
subset<-subset(data,GSI>0.6)

write.table(data.frame(table(subset$group)),file="Table.GSI.signature.for.each.tissue.txt",sep="\t",quote=F)

f3<-f2[which(gsi>0.3),]
colnames(f3)=colnames(file1)
colnames(f3)=colnames(f2)

histFfile=paste("Figure.Cluster.Dendrogram.plot.average",infile,"pdf",sep=".")
d <- dist(t(f3), method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()

# cluster and heatamp
newdata<-f2[which(gsi>0.6),]
pdf("Figure.supervised.histplot.analysis.hml.combat.quantile.pdf")
hist(newdata,breaks=40)
dev.off()

# yellow to blue
newdata[newdata<0]<-0
library("grDevices")
library("gplots")
pdf("Figure.supervised.heatmap.analysis.hml.combat.quantile.pdf")
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(newdata,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.65,labRow=NA)
dev.off()

# blue to yellow
newdata[newdata<0]<-0
library("grDevices")
library("gplots")
pdf("Figure.supervised.heatmap.analysis.hml.combat.quantile.pdf")
col=colorRampPalette(c("blue", "yellow"))(20) 
heatmap.2(newdata,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.65,labRow=NA)
dev.off()

# heatamp without cluster
data<-read.table(file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
pdf("GSI.hist.distribution.pdf")
hist(data$GSI,breaks=500)
dev.off()

# heat only have 3 high GSI regions
sum(table(subset(data,GSI>0.6)[,2]))
table(subset(data,GSI>0.6)[,2])
head(data)
tissue<-names(table(data[,2]))
tissue<-sort(c("Brain","Heart","muscle","Vessel","Spleen","Kidney","Ovary","Esophagus","Thymus","Lung","Liver","Pancreas","Stomach","Gastric","Intestine","Colon","Bladder"))
tissue<-names(table(colnames(f2)))[-c(3,7,8,21)]
  
choose<-c()
for(i in 1:length(tissue)){
  tmp<-subset(data,group==tissue[i])
  tmp<-tmp[order(tmp[,3],decreasing=T),]
  if(nrow(tmp)>5){
    choose<-c(choose,tmp[1:80,1])  
  }else{
    choose<-c(choose,tmp[,1])    
  }
  choose
}

tmp2<-f2[match(choose,rownames(f2)),]
tmp2<-tmp2[,unlist(lapply(tissue,function(x) grep(x,colnames(tmp2))))]
tmp2[tmp2<0]<-0.3
write.table(tmp2,file="high.gsi.tissue.matrix.txt",sep='\t',quote=F,col.names=NA,row.names=T)


library("grDevices")
library("gplots")
filename=paste("Figure.TSI",60,"pdf",sep=".")
pdf(filename)
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(tmp2,col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.8,labRow=NA)
dev.off()

tmp2[tmp2<0]<-0

# GWBS data analysis
setwd("/home/shg047/monod/nov")
load("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\nov\\GWBS.Data.RData")
infile="WGBS_methHap_load_matrix_16Oct2015.txt";
file1<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)

# preliminary clean

# sample name transform and normalized
colnames(file1) 
colnames(file1)<-gsub("_","-",colnames(file1))
colname2<-unlist(lapply(colnames(file1),function(x) unlist(strsplit(x,"[.]"))[1]))
colnames(file1)<-colname2
# be sure all the sample information has been stored in the following database
saminfo2<-read.table("/home/shg047/monod/phase2/newsaminfo.txt",head=T,sep="\t",as.is=T)
saminfo2<-saminfo2[match(colname2,saminfo2[,1]),]

# re-assign colnames
colnames(file1)<-saminfo2[,2]
colnames(file1)

# miss value imputation
library("impute")
f2<-RawNARemove(file1,missratio=0.3)
f2<-impute.knn(data.matrix(f2))$data
sum(is.na(file1))/(nrow(file1)*ncol(file1))
sum(is.na(f2))/(nrow(f2)*ncol(f2))

# check imptation 
ft<-file1[match(rownames(f2),rownames(file1)),]
f2[1:5,1:5] 
ft[1:5,1:5]

# quantile normalization
# library("preprocessCore")
# f3<-normalize.quantiles(f2,copy=TRUE)
# colnames(f3)<-colnames(f2)
# rownames(f3)<-rownames(f2)
# colnames(f3)

# batch effect (after clean and normalization)
# library("sva")
# batch=c(rep(1,2),rep(2,10),rep(3,10),rep(4,36),rep(5,3))
# colnames(f3)
# rownames(f3)

# check the boxplot distribution after combat
# save(f2,file="GWBS.Data.RData")
histFfile=paste("hist.combat",infile,"pdf",sep=".")
pdf(histFfile)
boxplot(f2,outline=F,horizontal=T,notch=F,las=1,cex.axis=0.65,col="blue")
dev.off()

# check the dendrogram after missing, quantile and combat
histFfile=paste("Cluster.Dendrogram.plot",infile,"pdf",sep=".")
f3 <- t(na.omit(f2)) # listwise deletion of missing
f3 <- scale(f3) # standardize variables
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
pdf(histFfile)
plot(fit,cex=0.5) # display dendogram
dev.off()


# check the boxplot distribution
histFfile=paste("hist",infile,"pdf",sep=".")
pdf(histFfile)
boxplot(f2,cex=0.5)
dev.off()

# check the cluster dendrogram
histFfile=paste("Cluster.Dendrogram.plot",infile,"pdf",sep=".")
f3 <- t(na.omit(f2)) # listwise deletion of missing
f3 <- scale(f3) # standardize variables
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
pdf(histFfile)
plot(fit,cex=0.5) # display dendogram
dev.off()

# check the cluster dendrogram with pvalue by bootstrap
library("pvclust")
fit <- pvclust(mydata, method.hclust="ward",method.dist="euclidean")
plot(fit)               # dendogram with p values
pvrect(fit, alpha=.95)  # add rectangles around groups highly supported by the data

# remove H1 and WBC, cancer
file1<-data.matrix(file1[,13:58])  

# miss value imputation
f2<-RawNARemove(file1,missratio=0.4)
file1<-impute.knn(f2)$data
dim(file1)

group=names(table(colnames(file1)))
index=colnames(file1)
gsi<-c()
gmaxgroup<-c()
for(i in 1:nrow(file1)){
  gsit<-0
  gmax<-names(which.max(tapply(as.numeric(file1[i,]),index,mean)))
  for(j in 1:length(group)){
    tmp<-(1-10^(mean(file1[i,][which(index==group[j])]))/10^(mean(file1[i,][which(index==gmax)])))/(length(group)-1)
    gsit<-gsit+tmp
  }
  gmaxgroup<-c(gmaxgroup,gmax)
  gsi<-c(gsi,gsit)
  print(c(gmax,gsit))
}
rlt=data.frame(region=rownames(file1),group=gmaxgroup,GSI=gsi)
write.table(rlt,file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",col.names=T,row.names=F,quote=F,sep="\t")

# each take top 5 tissue-specific methylation regions.
data<-read.table(file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
pdf("GSI.hist.distribution.pdf")
hist(data$GSI,breaks=200)
dev.off()
# heat only have 3 high GSI regions
sum(table(subset(data,GSI>0.6)[,2]))
table(subset(data,GSI>0.6)[,2])
head(data)
tissue<-names(table(data[,2]))
tissue<-sort(c("Brain","Heart","muscle","Vessel","Spleen","Kidney","Ovary","Esophagus","Thymus","Lung","Liver","Pancreas","Stomach","Gastric","Intestine","Colon","Bladder"))
choose<-c()
for(i in 1:length(tissue)){
  tmp<-subset(data,group==tissue[i])
  tmp<-tmp[order(tmp[,3],decreasing=T),]
  if(nrow(tmp)>5){
    choose<-c(choose,tmp[1:80,1])  
  }else{
    choose<-c(choose,tmp[,1])    
  }
  choose
}
tmp2<-file1[match(choose,rownames(file1)),]
tmp2<-tmp2[,unlist(lapply(tissue,function(x) grep(x,colnames(tmp2))))]
write.table(tmp2,file="high.gsi.tissue.matrix.txt",sep='\t',quote=F,col.names=NA,row.names=T)

library("grDevices")
library("gplots")
filename=paste("Figure.TSI",60,"pdf",sep=".")
pdf(filename)
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(tmp2,col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.8,labRow=NA)
dev.off()


# tissues signatures
names(table(data[,2]))
lung.signature<-subset(data,GSI>0.5 & group=="Lung")  # 0.52 for bspp
colon.signature<-subset(data,GSI>0.55 & group=="Colon")
pancrease.signature<-subset(data,GSI>0.68 & group=="Pancreas")
nrow(lung.signature)
nrow(colon.signature)
nrow(pancrease.signature)






################################################################################################
###########################################RRBS#################################################
################################################################################################

## cd /home/shg047/monod/hap/rrbs/bychrosome
## cp RRBS.mhl.matrix.txt /home/shg047/monod/dec/RRBS_methHap_load_matrix_Oct2015.txt
## cp RRBS.mhl.matrix.txt /home/shg047/monod/dec/RRBS_methaml_load_matrix_Oct2015.txt
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
hist(na_number,xlim=c(0,40))
dev.off()

# remove missing value and imputation

library("impute")
f2<-RawNARemove(file2,missratio=0.3)
f2<-impute.knn(data.matrix(f2))$data


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







