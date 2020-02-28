setwd("/home/sguo/drug/bamsample");
files<-list.files(pattern=".reads.normalize")
#source("fclust.R")

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

load("sample233.RData")
subs<-substr(files,1,12)
tmp<-match(unlist(subs),sample1[,1])
col<-data.matrix(data.frame(sample1[tmp,13]))-1


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
gene[is.na(gene)]<-0

#plot the gene expression of each gene
#each one has its own col
colu<-apply(gene,2,function(x) all(x==0))
if(any(colu)){
gene<-gene[,-which(colu==T)]
}

#each type(resist vs sensi) each line
pdfile=paste(list[i,4],".2.pdf",sep="")
pdf(pdfile)
sen<-apply(gene[which(col==1),],2,mean,trim = .2)
res<-apply(gene[which(col==2),],2,mean,trim = .2)
plot(sen,col="green",type="l")
lines(res,col="red")
legend("topright",legend=c("sen","res"),col=c("green","red"),lty=1,pch=1)
dev.off()

#each sample each line
#pdfile=paste(list[i,4],".pdf",sep="")
#pdf(pdfile)
#plot(gene[1,],type="l",col=col[1],ylim=c(0,100),main=list[i,4],ylab="Reads number of each site")
#for(m in 2:(length(files))){
#lines(gene[m,],type="l",col=col[m])
#}
#dev.off()
#}

}
