########################################################################################
###   Title: Heatamp plot for genome-wide BS-seq dataset from MoNOD Project
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 12/18/2015
########################################################################################

args = commandArgs(trailingOnly=TRUE)

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

library("impute")
################################################################################################
###########################################RRBS#################################################
################################################################################################
infile="WGBS_methHap_load_matrix_20Oct2015.txt";
file1<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)

infile="Heatmap.MHL.txt"
file2<-read.table(infile,head=T,sep="\t",row.names=1,as.is=T,check.names=F)
dim(file1)
dim(file2)

file1<-file1[match(rownames(file2),rownames(file1)),]
file1<-cbind(file1,file2)

# miss value detection and imputation
library("impute")
f2<-RawNARemove(file1,missratio=0.20)
f2<-impute.knn(data.matrix(f2))$data

library("preprocessCore")
# f2[,c(1:12,62:65)]<-normalize.quantiles(f2[,c(1:12,62:65)])
# f2[,13:61]<-normalize.quantiles(f2[,13:61])

library("sva")
# batch=c(rep(1,2),rep(2,10),rep(3,10),rep(4,36),rep(5,3),rep(6,4))
# f2<-ComBat(f2, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)
# save(f2,file="MHL.Trim.RData")

# load("MHL.Trim.RData")
# re-assign colnames
colnames(f2) 
colnames(f2)<-gsub("_","-",colnames(f2))
colname2<-unlist(lapply(colnames(f2),function(x) unlist(strsplit(x,"[.]"))[1]))
colname2
colnames(f2)<-colname2
# be sure all the sample information has been stored in the following database
saminfo2<-read.table("newsaminfo.txt",head=T,sep="\t",as.is=T)
saminfo2<-saminfo2[match(colname2,saminfo2[,1]),]
colnames(f2)<-saminfo2[,2]
colnames(f2)

saminfo3<-read.table("tissue2Layer.txt",head=T,sep="\t",as.is=T)
saminfo3<-saminfo3[match(colnames(f2),saminfo3[,1]),]
colnames(f2)<-saminfo3[,2]
colnames(f2)

file3 = f2
library(gplots)
library(RColorBrewer)
library("grDevices")
mydata <- file3
#colnames(mydata)[2]="CCT"
#GSI<-gsi(mydata)
#GSI<-GSI[order(GSI[,3],decreasing=T),]
# mydata<-mydata[match(names(sort(apply(mydata,1,sd),decreasing=T))[1:5000],rownames(mydata)),]
x<-apply(mydata,1,function(x) sd(x))
#y<-order(x,decreasing=T)[1:args[1]]
y<-order(x,decreasing=T)[1:args[1]]

mydata<-mydata[match(names(x[y]),rownames(mydata)),]
#mydata<-mydata[match(GSI[1:8000,1],rownames(mydata)),]
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")
# perform clustering on rows and columns
cl.row <- hclustfunc(distfunc(mydata))
cl.col <- hclustfunc(distfunc(t(mydata)))
# extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
gr.row <- cutree(cl.row, 6)
# gr.col <- cutree(cl.col, 5)
# require(RColorBrewer)
col1 <- brewer.pal(6, "Set1")     # the maximum of brewer.pal is 12
#col2 <- brewer.pal(5, "Pastel1")
#tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
tol21rainbow <- c("dodgerblue2","#E31A1C", # red
                  "green4",
                  "#6A3D9A", # purple
                  "#FF7F00", # orange
                  "black","gold1",
                  "skyblue2","#FB9A99", # lt pink
                  "palegreen2",
                  "#CAB2D6", # lt purple
                  "#FDBF6F", # lt orange
                  "gray70", "khaki2",
                  "maroon","orchid1","deeppink1","blue1","steelblue4",
                  "darkturquoise","green1","yellow4","yellow3",
                  "darkorange4","brown")
col2<-tol21rainbow[as.numeric(as.factor(cl.col$labels))]
col=colorRampPalette(c("blue", "yellow"))(20) 
require(gplots)    
filename=paste("heatmap_for_genomewide_mhl_topvar_sd_layer",args[1],"pdf",sep=".")
pdf(filename)
#par(mar=c(5,5,5,10))
heatmaprlt<-heatmap.2(as.matrix(mydata),hclustfun=hclustfunc, distfun=distfunc,
                      RowSideColors=col1[gr.row], 
                      ColSideColors=col2,
                      labRow=F,
                      trace="none",
                      col=col,
                      density.info="none")
#legend=unique(data.frame(col2,cl.col$labels))
#legend(x=0.95,y=0.8,legend=legend[,2],col=as.character(legend[,1]),pch=15,cex = 0.5)
dev.off()
save.image()
