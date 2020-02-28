library("grDevices")
library("gplots")
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
gsi<-function(data){
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
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=GSI)
  return(rlt)
}

TopGSIByCategory<-function(gsi,top=150){
  GSIRlt<-c()
  group<-names(table(gsi$group))
  rank<-c(rep(top,length(group)))
  for (i in 1:length(group)){
    subset=gsi[which(gsi$group==group[i]),]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
  }
  return(na.omit(GSIRlt))
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

rename<-function(data){
  Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T|h1",colnames(data))]
  colnames(Data)<-gsub("fileID_","",colnames(Data))
  colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
  colnames(Data)<-gsub("[.]","-",colnames(Data))
  colnames(Data)[grep("age|new|centenarian",colnames(Data))]<-"WBC"
  colnames(Data)[grep("X7.T|X6.T|SRX|CTT",colnames(Data))]<-"CT"
  colnames(Data)[grep("h1",colnames(Data))]<-"H1"
  colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
  return(Data)
}

frm<-function(data){
  # feature reduction (WGBS,RRBS and each caterogy missing<60,Plasma missing<50%,low variation removed)
  # features which are missing in whole reference samples will be omit since this will caused miss-classification (N<=4,data qualitye dependent)
  rm1<-which(apply(data[,grep("X7.P|X6.P|.6P|.7P|NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.5))
  rm2<-which(apply(data[,grep("X6.P|.6P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
  rm3<-which(apply(data[,grep("X7.P|.7P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
  rm4<-which(apply(data[,grep("NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
  rm5<-which(apply(data,1,function(x) sum(is.na(x))/length(x)>0.6))
  rm<-unique(c(rm1,rm2,rm3,rm4,rm5))
  return(rm)
}

# Dinh show me her data but I want to use my own dataset
# setwd("/media/Home_Raid1/dinh/WGBS_HaploInfo/scripts/Sharing/")
# infile1="wgbs.mhl.matrix";
# file1<-read.table(infile1,head=T,sep="\t",row.names=1,as.is=T,check.names=F)
# samplename<-saminfo[match(unlist(lapply(colnames(file1),function(x) gsub("fileID_","",x))),saminfo[,1]),2]
# colnames(file1)<-samplename
# setwd("/home/shg047/oasis/monod/hapinfo")
# scp /media/Home_Raid1/dinh/WGBS_HaploInfo/scripts/Sharing/wgbs.imf.matrix shg047@genome-miner.ucsd.edu:/media/NAS3_volume2/shg047/monod/hapinfo/
# Dinh's IMF format is better since she prepared MHB before CpG position

setwd("/media/Home_Raid1/shg047/work/monod/hapinfo")
library("impute")
library("preprocessCore")
library("sva")
saminfo<-read.table("N37Salk.saminfo",sep="\t")
data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
dataAMF<-read.table("AMF4.txt",head=T,row.names=1,sep="\t")
dataIMF<-read.table("/media/Home_Raid1/dinh/WGBS_HaploInfo/scripts/Sharing/wgbs.imf.matrix",head=T,row.names=1,sep="\t")
save(data,file="MHL4.RDATA")
save(dataAMF,file="AMF4.RDATA")
save(dataIMF,file="IMF4.RDATA")
load("MHL4.RDATA")
load("AMF4.RDATA")
load("IMF4.RDATA")

figure1A<-function(data,saminfo){
  data<-data[,grep("STL|N37|SRX|_h1|age|new|centenarian|CTT|HCT",colnames(data))]
  samname<-colnames(data)
  data<-rename(data)
  apply(data,2,function(x) sum(is.na(x)))
  data1<-RawNARemove(data,missratio=0.3)
  apply(data1,2,function(x) sum(is.na(x)))
  data2<-ColNARemove(data1,missratio=0.3)
  dim(data2)
  data3<-impute.knn(data.matrix(data2))$data
  colnames(data3)<-unlist(lapply(colnames(data3),function(x) unlist(strsplit(x,"[.]"))[1]))
  for(i in colnames(data3)){
  data3[,grep(i,colnames(data3))]<-normalize.quantiles(data.matrix(data3[,grep(i,colnames(data3))]))
  }
  #batch=c(rep(1,10),rep(2,36))
  #data4<-ComBat(data3, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)
  #data2[,1:10]<-normalize.quantiles(data.matrix(data3[,1:10]))
  #data2[,11:36]<-normalize.quantiles(data.matrix(data2[,11:36]))
  data4<-data3[order(apply(data3,1,sd),decreasing=T)[1:round(0.15*nrow(data3))],]
  png("Figure1C-NonSupervisor.png")
  col=colorRampPalette(c("blue", "yellow"))(20) 
  heatmap.2(data4,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.65,labRow=NA)
  dev.off()
  # non-supervisor
  GSIrlt<-gsi(data3)
  topgsi<-subset(GSIrlt,GSI>0.1)
  data4<-data3[match(topgsi[,1],rownames(data3)),]
  png("Figure1C-GSISupervisor.png")
  col=colorRampPalette(c("yellow", "blue"))(20) 
  heatmap.2(data4,col=col,trace="none",density.info="none",Colv=T,Rowv=T,key=T,keysize=1,cexCol=0.65,labRow=NA)
  dev.off()
  
  dataAMF<-dataAMF[,grep("STL|N37",colnames(dataAMF))]
  samname<-colnames(dataAMF)
  dataAMF<-rename(dataAMF)
  dataAMF2<-dataAMF[match(rownames(data4),rownames(dataAMF)),order(colnames(dataAMF))]
  dataAMF3<-impute.knn(data.matrix(dataAMF2))$data
  png("Figure1a.png")
  col=colorRampPalette(c("yellow", "blue"))(5) 
  heatmap.2(data.matrix(dataAMF3),col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.65,labRow=NA)
  dev.off()
}


figure1C<-function(data,saminfo){
  data<-data[,grep("STL|N37",colnames(data))]
  samname<-colnames(data)
  data<-rename(data)
  apply(data,2,function(x) sum(is.na(x)))
  data1<-RawNARemove(data,missratio=0.6)
  apply(data1,2,function(x) sum(is.na(x)))
  data2<-ColNARemove(data1,missratio=0.6)
  dim(data2)
  data3<-impute.knn(data.matrix(data2))$data
  colnames(data3)<-unlist(lapply(colnames(data3),function(x) unlist(strsplit(x,"[.]"))[1]))
  GSI<-gsi(data3)
  topgsi<-TopGSIByCategory(GSI,20)
  data4<-data3[match(topgsi[,1],rownames(data3)),]
  data4<-data4[,order(colnames(data4))]
  #batch=c(rep(1,10),rep(2,36))
  #data4<-ComBat(data3, batch, mod=NULL, par.prior = TRUE,prior.plots = FALSE)
  #data2[,1:10]<-normalize.quantiles(data.matrix(data3[,1:10]))
  #data2[,11:36]<-normalize.quantiles(data.matrix(data2[,11:36]))
  png("Figure1C-1.png")
  col=colorRampPalette(c("yellow", "blue"))(20) 
  heatmap.2(data4,col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.65,labRow=NA)
  dev.off()
  
  # Figure 1C-2
  dataAMF<-dataAMF[,grep("STL|N37",colnames(dataAMF))]
  samname<-colnames(dataAMF)
  dataAMF<-rename(dataAMF)
  dataAMF2<-dataAMF[match(rownames(data4),rownames(dataAMF)),order(colnames(dataAMF))]
  dataAMF3<-impute.knn(data.matrix(dataAMF2))$data
  png("Figure1C-2.png")
  col=colorRampPalette(c("yellow", "blue"))(5) 
  heatmap.2(data.matrix(dataAMF2),col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.65,labRow=NA)
  dev.off()
  
  # Figure 1C-3
  dataIMF<-dataIMF[,grep("STL|N37",colnames(dataIMF))]
  samname<-colnames(dataIMF)
  dataIMF<-rename(dataIMF)
  dataIMF2<-dataIMF[unlist(lapply(rownames(data4),function(x) grep(x,rownames(dataIMF)))),]
  dataIMF2<-impute.knn(data.matrix(dataIMF2))$data
  png("Figure1C-3.png")
  col=colorRampPalette(c("yellow", "blue"))(5) 
  heatmap.2(data.matrix(dataIMF2),col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.65,labRow=NA)
  dev.off()
}
