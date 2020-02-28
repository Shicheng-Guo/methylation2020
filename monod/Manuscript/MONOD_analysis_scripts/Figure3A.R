########################################################################################
###   Title: Heatamp plot for genome-wide BS-seq dataset from MoNOD Project
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 12/18/2015
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

gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,mean)))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(data[i,][which(index==group[j])]))/10^(mean(data[i,][which(index==gmax)])))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    gsi<-c(gsi,gsit)
    print(c(gmax,gsit))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi)
  return(rlt)
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

rename<-function(data){
  Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|h1",colnames(data))]
  colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
  colnames(Data)<-gsub("[.]","-",colnames(Data))
  colnames(Data)[grep("age|new|centenarian",colnames(Data))]<-"WBC"
  colnames(Data)[grep("X6.T|SRX381569|CTT",colnames(Data))]<-"LCT"
  colnames(Data)[grep("X7.T|SRX381716|SRX381719|SRX381722",colnames(Data))]<-"LCT"
  colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
  colnames(Data)[grep("X6.P|RRBS.6P",colnames(Data))]<-"CCP"
  colnames(Data)[grep("X7.P|RRBS.7P",colnames(Data))]<-"LCP"
  colnames(Data)[grep("h1",colnames(Data))]<-"H1"
  return(Data)
}

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}   

################################################################################################
###########################################RRBS#################################################
################################################################################################
setwd("/home/shg047/oasis/monod/hapinfo")
file1<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
saminfo<-read.table("/home/shg047/oasis/monod/saminfo/N37Salk.saminfo",sep="\t")
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
library("impute")
#f2<-RawNARemove(file1,missratio=0.20)
#f3<-ColNARemove(f2,missratio=0.60)
#f4<-impute.knn(data.matrix(f3))$data
file1<-file1[,grep("STL|N37|SRX|age|new|centenarian|CTT|HCT|h1",colnames(file1))]
f2<-RawNARemove(file1,missratio=0.60)
f2<-impute.knn(data.matrix(f2))$data
colnames(f2)

library("preprocessCore")
f3<-f2
#f3[,c(1:12,54:63)]<-normalize.quantiles(f2[,c(1:12,54:63)])
#f3[,c(13:52)]<-normalize.quantiles(f2[,c(13:52)])
#f3[,c(53,64:65)]<-normalize.quantiles(f2[,c(53,64:65)])

library("sva")
batch=c(rep(1,2),rep(2,10),rep(3,4),rep(4,36),rep(5,1),rep(1,10),rep(5,2))
f3<-ComBat(f3, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)

# re-assign colnames
f4<-rename(f3)
colnames(f4)<-unlist(lapply(colnames(f4),function(x) unlist(strsplit(x,"[.]"))[1]))

library(gplots)
library(RColorBrewer)
library("grDevices")
mydata <- f4
GSI<-gsi(mydata)
GSI<-GSI[order(GSI[,3],decreasing=T),]
# mydata<-mydata[match(names(sort(apply(mydata,1,sd),decreasing=T))[1:5000],rownames(mydata)),]


mydata<-mydata[match(GSI[1:5000,1],rownames(mydata)),]
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
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
col2<-tol21rainbow[as.numeric(as.factor(cl.col$labels))]
col=colorRampPalette(c("blue", "yellow"))(20) 
require(gplots)    
pdf("Figure3A.pdf")
heatmaprlt<-heatmap.2(as.matrix(mydata),hclustfun=hclustfunc, distfun=distfunc,
          RowSideColors=col1[gr.row], 
          ColSideColors=col2,
          labRow=F,
          trace="none",
          col=col,
          density.info="none")
legend=unique(data.frame(col2,cl.col$labels))
legend(x=0.85,y=0.8,legend=legend[,2],col=as.character(legend[,1]),pch=15,cex = 0.5)
dev.off()



