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


n <- 10000
x1  <- matrix(rnorm(n), ncol = 2)
x2  <- matrix(rnorm(n, mean = 3, sd = 1.5), ncol = 2)
x   <- rbind(x1, x2)

oldpar <- par(mfrow = c(2, 2))
smoothScatter(x, nrpoints = 0)
smoothScatter(x)


#########################################################################################

setwd("/home/sguo/monod/hap/wgbs/All_chromosomes_combined")
infile="wgbs.mhl.txt";
file1<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)
library("impute")
f2<-RawNARemove(file1,missratio=0.01)
f2<-impute.knn(data.matrix(f2))$data

sum(is.na(file1))/(nrow(file1)*ncol(file1))
sum(is.na(f2))/(nrow(f2)*ncol(f2))

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
# f3<-ComBat(f3, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)
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