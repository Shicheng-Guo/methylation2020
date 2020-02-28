ggbarplot<-function(data){
library("ggplot2")
newdata<-data.frame(value=as.numeric(data.matrix(data)),Group=rep(colnames(data),each=nrow(data)))
myData <- aggregate(newdata$value,by =list(type=newdata$Group),
                    FUN = function(x) c(mean = mean(x,na.rm=T),
                                        sd = sd(x,na.rm=T),
                                        sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                        min=min(x,na.rm=T),
                                        max=max(x,na.rm=T),
                                        median=median(x,na.rm=T),
                                        me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
)
myData <- do.call(data.frame, myData)
colnames(myData)=c("type","mean","sd","sem","min","max","median","me")
ggplot<-ggplot(myData, aes(x =type, y = mean)) +
  geom_bar(position = position_dodge(), stat="identity", fill="blue",width=0.8) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.0) +
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  coord_flip()+
  theme(axis.text=element_text(size=10),axis.title=element_text(),axis.text.y = element_text(hjust=0))
return(ggplot)
}

ggboxplot<-function(data){
library("ggplot2")
newdata<-data.frame(value=as.numeric(data.matrix(data)),Group=rep(colnames(data),each=nrow(data)))
ggplot<-ggplot(newdata, aes(factor(Group),value)) +
  geom_boxplot(aes(fill = factor(Group)),outlier.shape=NA)+
  coord_flip()+
  theme(axis.text=element_text(size=10),axis.title=element_text(),axis.text.y = element_text(hjust=0),legend.position="none")
return(ggplot)
}

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    dat<-data[,-NaCOL]
  }else{
    dat<-data;
  }
  dat
}   

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

ttestFunction<-function(data,x1,x2){
  data=data+matrix(rnorm(length(data),0.000001,0.000001),nrow(data),ncol(data))
  output<-matrix(NA,dim(data)[1],4)
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),length(na.omit(data[i,x1]))<2,length(na.omit(data[i,x2]))<2,all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]), na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-as.numeric((mean(as.numeric(data[i,x1]),na.rm=T)-mean(as.numeric(data[i,x2]),na.rm=T)))
      output[i,3]<-mean(as.numeric(data[i,x1]),na.rm=T)
      output[i,4]<-mean(as.numeric(data[i,x2]),na.rm=T)
      # print(i)
    }
  }
  
  rownames(output)=rownames(data)
  P.Adj<-p.adjust(output[,1],method="fdr")
  out<-data.frame(output[,1],P.Adj,output[,2:4])
  out<-na.omit(out)
  colnames(out)=c("Pvalue","FDR","Delta","MG1","MG2")
  return(out)
}


setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
load("MONOD-Apr6.MHL.RData")

## Miss Value Statistic
sum(is.na(data))/(nrow(data)*ncol(data))
na<-apply(data,2,function(x) sum(is.na(x)))
sort(na)
barplot(sort(na),horiz=T,las=1,cex.names=0.25,col="blue")
write.table(sort(na),file="missing.value.table.txt",sep="\t")

## re-group/collect the samples #colon cancer
Group<-colnames(data)
newdata<-data
colnames(newdata)<-Group
newdata<-data.frame(newdata,check.names=F)

gsi<-function(data){
  # NA robust
  group=names(table(colnames(data)))
  index=colnames(data)
  GSI<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    GSI<-c(GSI,gsit)
    #print(paste(gmax,gsit,sep=" "))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=GSI)
  return(rlt)
}

PCAPlot<-function(data,pheno,output,multifigure=F){
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
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",pt.lwd=1.5)
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",pt.lwd=1.5)
  dev.off()
}
PCAPlot(data,pheno,output,multifigure=F)
load("MONOD-Apr6.MHL.RData")
Group<-colnames(data)
YY<-grep("STL|N37-|methylC-|adenocarcinoma_lung|tumor_lung|Colon_Tumor_Primary|CTT-|metastasis_colon|tumor_colon",Group)
Newdata<-data[,YY]
colnames(Newdata)
write.table(colnames(Newdata),file="PCA.Dataset.Colnames-tmp.txt",sep="\t")
# regroup of PCA dataset
saminfo<-read.table("PCA.Dataset.Colnames.txt",sep="\t",as.is=T)
Newdata<-Newdata[,colnames(Newdata) %in% saminfo[,1]]
pheno<-saminfo[match(colnames(Newdata),saminfo[,1]),2]
## Signficant Detection without remove missing values
Newdata<-RawNARemove(Newdata,0.3)
dim(Newdata)
Newdata[1:4,1:4]
library("impute")
Newdata<-impute.knn(data.matrix(Newdata))$data
# PCA based on total Dataset
data=t(Newdata) # row individual, column is gene
pheno=pheno
output="PCA1"
multifigure=F
PCAPlot(data,pheno,output,multifigure=F)
dim(data)
# PCA based on high variable data
NewdataGSI<-Newdata
colnames(NewdataGSI)=pheno
GSI<-gsi(NewdataGSI)
dim(GSI)
subset<-subset(GSI,GSI>0.3)  # to select high variable regions.
dim(subset)
head(subset)
NewData<-Newdata[match(subset[,1],rownames(Newdata)),]
dim(NewData)
NewData[1:4,1:4]
data=t(NewData) # row individual, column is gene
pheno=pheno
multifigure=F
PCAPlot(data,pheno,output="PCA2",multifigure=F)

