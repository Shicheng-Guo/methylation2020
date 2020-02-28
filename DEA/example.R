setwd("/home/sguo/drug");
files<-list.files(pattern=".reads.normalize")
source("fclust.R")

#setup list
gene.list<-read.table("drug_mnotu_feature.bed",head=F,sep="\t",as.is=F)
column.list<-c(0,cumsum(gene.list[,3]-gene.list[,2]+1)) # the first 6 column are fixed for ped file
len<-length(column.list)
list<-array()
for (i in 1:(len-1)){
  tmp<-paste(column.list[i]+1,column.list[i+1],i,sep="\t")
  list[i]<-tmp
}
xlist<-cbind(gene.list,list)
write.table(xlist,file="column.list.txt",col.names=F,row.names=F,sep="\t",quote=F)
list<-read.table("column.list.txt",sep="\t")
head(list)
#
i<-j<-numeric()
fpca.rlt<-list()
fpca.rnaseq.rlt<-matrix(NA,length(files))

for (i in 1:dim(list)[1]){
print (i);
skip<-list[i,7]-1
nrows<-list[i,3]-list[i,2]+1
ncols<-length(files)
gene<-matrix(rep(NA,nrows*ncols),nrow=nrows,ncol=ncols)
for (j in 1:ncols){
gene[,j]<-unlist(read.table(files[j],skip=skip,nrows=nrows))
}
colnames(gene)=files;
gene<-t(gene)
gene <- na.omit(gene) # listwise deletion of missing
#gene <- scale(gene) # standardize variables 
gene[is.na(gene)]<-0

#plot the gene expression of each gene
#col<-apply(gene,2,function(x) all(x==0))
#if(any(col)){
#gene<-gene[,-which(col==T)]
#}
#pdfile=paste(list[i,4],".pdf",sep="")
#pdf(pdfile)
#plot(gene[1,],type="l",col=1,ylim=c(0,50),main=list[i,4],ylab="Reads number of each site")
#for(m in 2:100){
#lines(gene[m,],type="l",col=m)
#}
#dev.off()
##

if(dim(gene)[2]>500){
filename=paste(list[i,4],i,sep="_")
fpca.rlt<- try( fpca.genotype(gene,pos=NULL,nbasis=21)  )#combine principle component
fpca.rnaseq.rlt<-cbind(fpca.rnaseq.rlt,fpca.rlt$scores);
print(dim(fpca.rnaseq.rlt))
save(fpca.rnaseq.rlt,file="fpca.result.noscale.RData")
print (i);
#write.table(gene,file=filename,col.names=NA,row.names=T,sep="\t")
}
}

x<-fpca.rnaseq.rlt[,-1]
y<-fhclust(x,method="euclidean",cut=2,bootstrap=FALSE) 
y 
write.table(x,file="drug.pca.noscale.txt",sep="\t",col.names=F,row.names=T,quote=F)
write.table(y,file="cluster.txt",sep="\t",col.names=F,row.names=T,quote=F)



