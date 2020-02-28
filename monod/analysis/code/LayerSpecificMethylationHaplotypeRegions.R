
########################################################################################
###   Title: layer specific methylation haplotype
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 11/9/2015
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

###################################################################################################################



data<-read.table(file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
# heat only have 3 high GSI regions


file1<-read.table("RRBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file1)
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
remove=c("6-T-3","6-T-4","7-T-2",paste("NC-P-",19:24,sep=""),"PC-P-10","6-P-6",paste("PC-P-",c(2,3,6,9),sep=""))
file1<-file1[,-match(remove,samplename2)]
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
new<-read.table("saminfo.txt",sep="\t",as.is=T)
cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
matrix=file1
samplename2<-gsub("6-P","CC-P",samplename2)
samplename2<-gsub("7-P","LC-P",samplename2)
samplename2<-gsub("6-T","CC-T",samplename2)
samplename2<-gsub("7-T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
colnames(matrix)=samplename2



lung.signature<-subset(data,GSI>0.689 & group=="Lung")   # number=10
colon.signature<-subset(data,GSI>0.7709 & group=="Colon") # number=10
pancrease.signature<-subset(data,GSI>0.8474 & group=="Pancreas") # number=10

lung.signature<-subset(data,GSI>0.66 & group=="Lung")
colon.signature<-subset(data,GSI>0.75 & group=="Colon")
pancrease.signature<-subset(data,GSI>0.79 & group=="Pancreas")

lung.signature<-subset(data,GSI>0.66 & group=="Lung")
colon.signature<-subset(data,GSI>0.75 & group=="Colon")
pancrease.signature<-subset(data,GSI>0.83 & group=="Pancreas")

nrow(lung.signature) 
nrow(colon.signature)
nrow(pancrease.signature)

bed1<-cor2bed(rownames(matrix))
bed2<-cor2bed(lung.signature[,1])
bed3<-cor2bed(colon.signature[,1])
bed4<-cor2bed(pancrease.signature[,1])

rlt.lung<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed2,opt.string="-wa -u")
rlt.colon<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed3,opt.string="-wa -u")
rlt.pancrease<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed4,opt.string="-wa -u")

cor.lung<-bed2cor(rlt.lung)
cor.colon<-bed2cor(rlt.colon)
cor.pancrease<-bed2cor(rlt.pancrease)

length(cor.lung)
length(cor.colon)
length(cor.pancrease)

data.prediction.lung<-file1[match(cor.lung,rownames(file1)),c(grep("6-P",colnames(file1)))]
data.prediction.colon<-file1[match(cor.lung,rownames(file1)),c(grep("7-P",colnames(file1)))]
data.prediction.pancreatic<-file1[match(cor.lung,rownames(file1)),c(grep("PC-P",colnames(file1)))]

t1<-apply(data.prediction.lung,2,function(x) sum(x>0)/length(x))
t2<-apply(data.prediction.colon,2,function(x) sum(x>0)/length(x))
t3<-apply(data.prediction.pancreatic,2,function(x) sum(x>0)/length(x))


#
choose.prediction<-unique(c(cor.lung,cor.colon,cor.pancrease))
data.prediction.data<-file1[match(choose.prediction,rownames(file1)),c(grep("6-P",colnames(file1)),grep("7-P",colnames(file1)),grep("PC-P",colnames(file1)))]
t<-apply(data.prediction.data,2,function(x) { u=x>0;print(c(sum(unlist(u)[1:10]),sum(unlist(u)[11:20]),sum(unlist(u)[21:30])))})








