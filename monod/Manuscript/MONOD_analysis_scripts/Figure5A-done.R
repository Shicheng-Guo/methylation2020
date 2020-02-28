
## 2017-05-12

source("http://www.bioconductor.org/biocLite.R")
biocLite("impute")  
#install.packages("gplots")
#install.packages("RColorBrewer")
#install.packages("grDevices")
library("gplots")
library("RColorBrewer")
library("grDevices")
library("impute")
setwd("/home/shg047/oasis/monod/hapinfo")
load("MHL4.RData")
data1<-data[,c(grep(".6P|X6.P",colnames(data)),grep("X6.T|CTT|SRX381569",colnames(data)),grep("Colon",colnames(data)),grep("NC.P|WB|age|born|centenarian",colnames(data)))]
x1<-c(grep(".6P|X6.P",colnames(data1)))
x2<-c(grep("X6.T|CTT|SRX381569|colon|Colon|HCT116|SG-",colnames(data1)),grep("Colon",colnames(data1)))
x3<-c(grep("NC.P|WB|age|born|centenarian",colnames(data1)))
length(x1)
length(x2)
length(x3)
# pre-select and then random-forest
data2<-data1[apply(data1,1,function(x) sum(x[x3]<0.1,na.rm=T)>0.75*length(na.omit(x[x3]))),]
tmp<-data2[,x3]
#tmp[tmp>0.2]<-0.01
data2[,x3]<-tmp
data2<-data2[which(apply(data2,1,function(x) sum(x[x2]>0.2)>=0.5*length(na.omit(x[x2])))),]  
data2<-data2[which(apply(data2,1,function(x) sum(x[x1]>0.01)>=0.2*length(na.omit(x[x1])))),]  
dim(data2)
mydata<-data2
dim(mydata)
f2<-RawNARemove(mydata,missratio=0.3)
mydata<-impute.knn(data.matrix(f2))$data
hclustfunc <- function(x) hclust(x, method="ward.D")
distfunc <- function(x) dist(x, method="euclidean")
cl.row <- hclustfunc(distfunc(mydata))
cl.col <- hclustfunc(distfunc(t(mydata)))
gr.row <- cutree(cl.row, 6)
col1 <- brewer.pal(6, "Set1") 
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
col2<-tol21rainbow[as.numeric(as.factor(cl.col$labels))]
filename=paste("Figure 5A.Colon.tissue.cancer.signature","pdf",sep=".")
pdf(filename)
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(mydata,col=col,trace="none",density.info="none",Colv=F,Rowv=T,key=T,keysize=1,cexCol=0.4,labRow=T,na.rm=TRUE)
dev.off()
getwd()
write.table(mydata,file="colon.tissue-Figure5A.signiture.txt",col.names=NA,row.names=T,sep="\t")
