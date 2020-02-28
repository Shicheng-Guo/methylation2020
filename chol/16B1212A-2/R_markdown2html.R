library("knitr")
library("rmarkdown")

data<-read.table("he2019.fam")
data[grep("RA",data$V2),6]<-1
data[-grep("RA",data$V2),6]<-0
write.table(data,file="he2019.fam",col.names = F,row.names = F,quote = F)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/SLE/BCR/vdj")
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ASCs/master/extdata/IFI.txt",head=T,sep="\t")
head(data)
input<-data[,8:ncol(data)]
head(input)
rownames(input)<-data[,2]
idx<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[.]"))[1]))
par(mfrow=c(3,4),mar=c(2,2,3,1))
for(i in 1:16){
  RNA<-log(as.numeric(input[i,]),2)
  boxplot(RNA~idx,col=2:3,main=rownames(input)[i])
}

data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ASCs/master/extdata/vdjtools.basicstats.txt",head=T,sep="\t")
head(data)
idx1<-unlist(lapply(data[,1],function(x) substr(x,1,1)))
idx2<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[_]"))[2]))
idx<-paste(idx1,idx2,sep="_")
dim(data)
par(mfrow=c(3,4),mar=c(2,2,3,1))
for(i in 3:ncol(data)){
  boxplot(data[,i]~idx,col=2:5,main=colnames(data)[i])
}

data<-read.table("vdjtools.segments.wt.J.txt",head=T,sep="\t")
head(data)
idx1<-unlist(lapply(data[,1],function(x) substr(x,1,1)))
idx2<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[_]"))[2]))
idx<-paste(idx1,idx2,sep="_")
dim(data)
par(mfrow=c(2,4),mar=c(2,2,3,1))
for(i in 3:ncol(data)){
  boxplot(data[,i]~idx,col=2:5,main=colnames(data)[i])
}

data<-read.table("vdjtools.segments.wt.V.txt",head=T,sep="\t")
head(data)
dim(data)
idx1<-unlist(lapply(data[,1],function(x) substr(x,1,1)))
idx2<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[_]"))[2]))
idx<-paste(idx1,idx2,sep="_")
dim(data)
par(mfrow=c(4,4),mar=c(2,2,3,1))
for(i in 3:ncol(data)){
  if(sd(data[,i])>0.0055){
    boxplot(data[,i]~idx,col=2:5,main=colnames(data)[i],ylim=c(0,0.12))
  }
}

