setwd("/home/sguo/Dropbox/Project/methylation/monod/phase2/")
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\phase2")
setwd("/home/sguo/Dropbox/Project/methylation/monod/phase2")
new<-read.table("saminfo.txt",sep="\t",as.is=T)
head(new)

file1<-read.table("RRBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file1)
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
remove=c("6-T-3","6-T-4","7-T-2",paste("NC-P-",19:24,sep=""),"PC-P-10","6-P-6",paste("PC-P-",c(2,3,6,9),sep=""))
file1<-file1[,-match(remove,samplename2)]
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
matrix=file1
samplename2<-gsub("6-P","CC-P",samplename2)
samplename2<-gsub("7-P","LC-P",samplename2)
samplename2<-gsub("6-T","CC-T",samplename2)
samplename2<-gsub("7-T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
colnames(matrix)=samplename2
d <- dist(t(matrix)) # distance matrix
fit <- hclust(d, method="complete")         # distance matrix
pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.7,hang=-1,xlab="",,ylab="",lwd=2.5,main="",cex.axis=0.7)

dev.off()
par(mar=c(3,5,3,3))
boxplot(matrix[,1:ncol(matrix)],outline=F,horizontal=T,notch=T,las=1,cex.axis=0.65)

par(mar=c(3,5,1,1))
newmatrix<-matrix[,order(apply(matrix,2,function(x) quantile(x)[4]),decreasing=T)]
boxplot(newmatrix[,1:ncol(newmatrix)],outline=F,horizontal=T,las=1,cex.axis=0.65)

file1<-read.table("WGBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file1)
colnames(file1)<-gsub("_","-",colnames(file1))
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
matrix=file1
samplename2<-gsub("6-P","CC-P",samplename2)
samplename2<-gsub("7-P","LC-P",samplename2)
samplename2<-gsub("6-T","CC-T",samplename2)
samplename2<-gsub("7-T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
colnames(matrix)=samplename2
matrix<-matrix[,-c(11,12)]
d <- dist(t(matrix)) # distance matrix
fit <- hclust(d, method="complete")         # distance matrix

pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.7,hang=-1,xlab="",,ylab="",lwd=2.5,main="",cex.axis=0.7)
dev.off()

dim(matrix)
pdf("wgbs.boxplot.pdf")
par(mar=c(2,8,2,2))
par(mai=c(2,15,2,2))

boxplot(matrix[,1:ncol(matrix)],outline=F,horizontal=T,notch=T,las=1,cex.axis=0.65)
dev.off()


phen=c(rep("Cancer",2),rep("H1",10),rep("N37-Normal",10),rep("Salk-Normal",36))
PCAPlot(data=t(matrix),pheno=phen,output="WGBS.PCA.Phase2",multifigure=T)
  
colnames(matrix)
phen=c(rep("Cancer",2),rep("H1",10),rep("N37-Normal",10),rep("Salk-Normal",36),"WB-CEN","WB-MID","WB-NEW")
data=t(matrix)
pheno=phen
output="WGBS.PCA.Phase2"
multifigure=T

  
  

file3<-read.table("WGBS_SeqCap_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file3)
samplename1=sapply(strsplit(colnames(file3),"_"),function(x) unlist(x)[1])



pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.7,hang=-1,xlab="",,ylab="",lwd=3,main="",cex.axis=0.7)
dev.off()



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
    points(scores$PC1[i]+rnorm(1,1,1),scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=1,lwd=2)
  }
  legend("topleft",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",cex=0.9)
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("topleft",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",cex=0.7)
  dev.off()
}

