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
setwd("/home/shg047/monod/wgbs")
saminfo<-read.table("/home/shg047/monod/saminfo/N37Salk.saminfo",head=T,sep="\t")
wgbsMf<-read.table("WGBS.N37.Normal.CpG.methfreq.matrix.txt",head=T,sep="\t",check.names=F,as.is=T,row.names=1)
highgsiBed<-read.table("/home/shg047/monod/dec/hgsi.wgbs.aml.bed",head=F,sep="\t")

library("impute")
f2<-RawNARemove(newWgbsMf,missratio=0.3)
f2<-impute.knn(data.matrix(f2))$data
colnames(f2)<-saminfo[match(colnames(f2),saminfo[,1]),2]
tissue<-names(table(colnames(f2)))
f2<-f2[,unlist(lapply(tissue,function(x) grep(x,colnames(f2))))]
write.table(f2,file="WGBS.N37.Normal.CpG.methfreq.matrix.txt",col.names=NA,row.names=T,quote=F,sep="\t")

bed1<-highgsiBed
bed2<-cor2bed(rownames(wgbs_mf))
bed<-Rbedtools(functionstring="intersectBed",bed1,bed2,opt.string="-wo")
newWgbsMf<-wgbsMf[match(bed2cor(bed[,4:6]),rownames(wgbsMf)),]
colnames(newWgbsMf)<-colnames(wgbs_mf)
newWgbsMf[1:5,1:5]

newdata<-newWgbsMf
newdata[newdata<0]<-0
library("grDevices")
library("gplots")
pdf("Figure.supervised.aml.single.cpg.heatmap.analysis.combat.quantile.pdf")
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(data.matrix(newdata),col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.65,labRow=NA)
dev.off()

bedtools intersect -wb -a /home/shg047/monod/dec/hgsi.wgbs.aml.bed -b /home/shg047/monod/wgbs/WGBS.N37.Normal.CpG.methfreq.matrix.txt > WGBS.N37.Normal.CpG.methfreq.matrix.intersect.bed