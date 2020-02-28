#######################################################################################################################
###   Title : Heatmap plot based on raw methylation signals which shared with high GSI MHL regions
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Time :  Sep/23/2015 
###   New: Extract the methylation signals with MethylFreq2Matrix.pl
#######################################################################################################################

data<-read.table("MethylMatrix.Filelist.txt.freq",row.names=1,head=T,sep="\t",as.is=T,check.names=F)
s1<-grep("KZ",colnames(data))
s2<-grep("STL",colnames(data))
length(s1)
length(s2)
newdata=data[,c(s1,s2)]
dim(na.omit(newdata))
newdata<-RawNARemove(newdata)
library("impute")
newdata<-impute.knn(data.matrix(newdata))$data
newdata[1:3,1:3]
saminfo<-read.table("saminfo.txt",head=F,sep="\t",as.is=T)
colnames(newdata)=saminfo[match(colnames(newdata),saminfo[,1]),2]
newdata<-newdata[,order(colnames(newdata))]

save(newdata,file="Raw.methy.GSI.RData")

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("grDevices")
# biocLite("gplots")

library("grDevices")
library("gplots")
pdf("Figure-20-raw-signal.pdf")
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(newdata,col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=NA,keysize=1,cexCol=0.8,labRow=NA)
dev.off()
newdata<-newdata[match(corsort(rownames(newdata)),rownames(newdata)),]


corsort<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) x))
  bed<-matrix(a,ncol=2,byrow=T)
  bed<-bed[order(bed[,1],as.numeric(bed[,2])),]
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],sep="")})
  return(cor)
}



