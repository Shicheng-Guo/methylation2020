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
for(i in 1:nrow(fn)){
  gsit<-0
  gmax<-names(which.max(tapply(as.numeric(fn[i,]),index,mean)))
  for(j in 1:length(group)){
    tmp<-(1-10^(mean(fn[i,][which(index==group[j])]))/10^(mean(fn[i,][which(index==gmax)])))/(length(group)-1)
    gsit<-gsit+tmp
  }
  gmaxgroup<-c(gmaxgroup,gmax)
  gsi<-c(gsi,gsit)
  print(c(gmax,gsit))
}
rlt=data.frame(region=rownames(fn),group=gmaxgroup,GSI=gsi)
write.table(rlt,file="Table.GSI.layer.mhl.WGBS.Remove.H1.WBC.rlt.txt",col.names=T,row.names=F,quote=F,sep="\t")

# heatamp without cluster
data<-read.table(file="Table.GSI.layer.mhl.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
pdf("Figure.layer.mhl.GSI.hist.distribution.pdf")
hist(data$GSI,breaks=200,xlim=c(0,1),col="red",cex=1.4)
dev.off()

subset<-subset(data,GSI>0.6)
newf2<-fn[match(subset[,1],rownames(fn)),]

newdata<-newf2
newdata<-newdata[,-c(11,12)]
newdata[newdata<0]<-0
library("grDevices")
library("gplots")
pdf("Figure.supervised.layer.mhl.single.cpg.heatmap.analysis.combat.quantile.pdf")
col=colorRampPalette(c("yellow", "blue"))(20) 
rlt<-heatmap.2(data.matrix(newdata),col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.75,labRow=NA)
dev.off()


write.table(newdata,file="Supplementary.layer.specific.mhl.matrix.txt",col.names=NA,row.names=T,quote=F,sep="\t")

newdata<-read.table("Supplementary.layer.specific.mhl.matrix.txt",head=T,sep="\t",row.names=1)
col=colorRampPalette(c("yellow", "blue"))(20) 
rlt<-heatmap.2(data.matrix(newdata),col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.65,labRow=NA)

hc<-hclust(dist(newdata))
clust<-cutree(hc,k=3)
k1<-names(clust[clust==1])
k2<-names(clust[clust==2])
k3<-names(clust[clust==3])

bed1<-cor2bed(k1)
bed2<-cor2bed(k2)
bed3<-cor2bed(k3)

dim(bed1)
dim(bed2)
dim(bed3)

write.table(bed1,file="ec.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(bed2,file="en.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(bed3,file="me.bed",sep="\t",quote=F,col.names=F,row.names=F)



## loop for layer specific
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





