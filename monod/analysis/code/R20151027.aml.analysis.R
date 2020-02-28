########################################################################################
###   Title: Average methylation level dataset analysis
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



###################################################################################################################

setwd("/home/shg047/monod/dec")
infile="WGBS_aveMeth_load_matrix_20Oct2015.txt";
file1<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)

# miss value detection and imputation
library("impute")
f2<-RawNARemove(file1,missratio=0.3)
f2<-impute.knn(data.matrix(f2))$data

library("preprocessCore")
f2.t1<-normalize.quantiles(f2[,13:58])
library("sva")
batch=c(rep(1,10),rep(2,36))
f2.t2<-ComBat(f2.t1, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)
f2[,13:58]<-f2.t2

# re-plot the cluster by sequencing sample Tag
colnames(f2)<-colnames(file1)
f3 <- t(na.omit(f2)) # listwise deletion of missing
histFfile=paste("Figure.Cluster.aml.Dendrogram.combat.seqTag.ward.D.plot",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()
histFfile=paste("Figure.Cluster.aml.Dendrogram.combat.seqTag.plot.average",infile,"pdf",sep=".")
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
# check the dendrogram after missing, quantile
f3 <- t(na.omit(f2)) # listwise deletion of missing
# f3 <- scale(f3) # standardize variables
sum(is.na(file1))/(nrow(file1)*ncol(file1))
sum(is.na(f2))/(nrow(f2)*ncol(f2))
histFfile=paste("Figure.aml.hist.combat",infile,"pdf",sep=".")
pdf(histFfile)
boxplot(f2,outline=F,horizontal=T,notch=F,las=1,cex.axis=0.65,col="blue")
dev.off()
histFfile=paste("Figure.Cluster.aml.tissuetag.Dendrogram.combat.ward.D.plot",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()
histFfile=paste("Figure.Cluster.aml.tissuetag.Dendrogram.combat.plot.average",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()


colnames(f2)<-colnames(file1)
# check the dendrogram after missing, quantile
f3 <- t(na.omit(f2)) # listwise deletion of missing
# f3 <- scale(f3) # standardize variables
sum(is.na(file1))/(nrow(file1)*ncol(file1))
sum(is.na(f2))/(nrow(f2)*ncol(f2))
histFfile=paste("Figure.aml.hist.combat",infile,"pdf",sep=".")
pdf(histFfile)
boxplot(f2,outline=F,horizontal=T,notch=F,las=1,cex.axis=0.65,col="blue")
dev.off()
histFfile=paste("Figure.Cluster.aml.seqtag.Dendrogram.combat.ward.D.plot",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()
histFfile=paste("Figure.Cluster.aml.seqtag.Dendrogram.combat.plot.average",infile,"pdf",sep=".")
d <- dist(f3, method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()


# tissue specific hml
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
write.table(rlt,file="Table.GSI.aml.WGBS.Remove.H1.WBC.rlt.txt",col.names=T,row.names=F,quote=F,sep="\t")

# heatamp without cluster
data<-read.table(file="Table.GSI.aml.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
pdf("Figure.aml.GSI.hist.distribution.pdf")
hist(data$GSI,breaks=200,xlim=c(0,1),col="pink")
dev.off()


f3<-f2[which(gsi>0.75),]
colnames(f3)=colnames(file1)
colnames(f3)=colnames(f2)

histFfile=paste("Figure.Cluster.aml.Dendrogram.plot.average",infile,"pdf",sep=".")
d <- dist(t(f3), method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()

# cluster and heatamp
newdata<-f2[which(gsi>0.85),]
pdf("Figure.supervised.aml.histplot.analysis.combat.quantile.pdf")
hist(newdata,breaks=40)
dev.off()

# yellow to blue
newdata[newdata<0]<-0
library("grDevices")
library("gplots")
pdf("Figure.supervised.aml.heatmap.analysis.combat.quantile.pdf")
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(newdata,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.65,labRow=NA)
dev.off()

# blue to yellow
newdata[newdata<0]<-0
library("grDevices")
library("gplots")
pdf("Figure.supervised.aml.heatmap.analysis.combat.quantile.pdf")
col=colorRampPalette(c("blue", "yellow"))(20) 
heatmap.2(newdata,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.65,labRow=NA)
dev.off()

# heatamp without cluster
data<-read.table(file="Table.GSI.aml.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
pdf("GSI.hist.distribution.pdf")
hist(data$GSI,breaks=500)
dev.off()



