#######################################################################################################################
###   Title : Heatmap plot based on raw methylation signals which shared with high GSI MHL regions
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Time :  Sep/23/2015 
###   New: Extract the methylation signals with MethylFreq2Matrix.pl
###   Prerequisite: achieve bed file (target region)
###   Prerequisite: MethylFreq files for all the sample
###   Directory: 512 server: /home/sguo/monod/methyFreq/
###   09/18/2015: MLH data updated and bioinformatic analysis again
#######################################################################################################################

library("Biostrings")
library("stringr")
library("randomForest")
library("impute")
library("rpart")
library("e1071")
library("biclust")

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}
pipeline<-function(lab1=lab1,file1=file1,cor1=cor1,dataset="RRBS"){
  rlt<-list()
  u<-which(lab1=="Cancer")
  v<-which(lab1=="Normal")
  rlt$cancer=u
  rlt$normal=v
  rlt$rawdata=file1
  keep<-c()
  for(i in 1:nrow(file1)){
    p1<-sum(file1[i,u]>0)
    p2<-sum(file1[i,v]==0)
    print(i)
    if(p2==length(v) & p1>0){
      tmp<-c(i,p1)
      keep<-rbind(keep,tmp)
    }  
  }
  rlt$totalnum<-nrow(file1)
  rlt$hyperfreq<-keep
  rlt$hypernum<-nrow(keep)
  
  # hypermethylation frequency in cancer samples for the loci which is totally un-methylation in normal samples. 
  Ofile1<-paste(dataset,"hypermethylation.freq.hist","jpeg",sep=".")
  jpeg(Ofile1)
  hist(keep[,2]/length(u),main="",xlab="Hypermethylation Frequency in Cancer",lwd=1.5,cex.lab=1.5,cex.axis=1.5)
  dev.off()
  
  # function used to get best mtry and mtree
  data<-as.matrix(file1[keep[,1],])
  bestmtry<-tuneRF(x=t(data), y=as.factor(lab1), ntreeTry=1000, stepFactor=2, improve=0.05,trace=TRUE, plot=TRUE, doBest=T)
  best.mtry=bestmtry$mtry
  senspeacc<-function(x){
    acc=(x$confusion[1,1]+x$confusion[2,2])/sum(x$confusion)
    sen=(x$confusion[1,1])/(x$confusion[1,1]+x$confusion[1,2])
    spe=(x$confusion[2,2])/(x$confusion[2,1]+x$confusion[2,2])
    c(sen,spe,acc)
  }
  err.rate<-c()
  senspacc<-c()
  for(mtrees in seq(250,5000,by=250)){
    rf.model1<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=mtrees)
    tmp1<-mean(rf.model1$err.rate[,1])
    tmp2<-senspeacc(rf.model1)
    err.rate<-c(err.rate,tmp1)
    senspacc<-c(senspacc,tmp2)
  }
  bestmtrees<-seq(250,5000,by=250)[which.min(err.rate)]
  
  # do the random forest based on above best mtry and mtree
  # informative loci which is useful in Random Forest model based on methylated regions in cancer tissues. 
  rf.model<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=bestmtrees,importance=FALSE)
  importance<-data.frame(ord=rownames(rf.model$importance),importance=rf.model$importance)
  num<-sum(importance[,2]>0)
  rlt$acc1<-senspeacc(rf.model)
  
  
  newdata1<-data[order(importance[,2],decreasing=T)[1:num],]
  keep1<-c()
  for(i in 1:nrow(newdata1)){
    p1<-sum(newdata1[i,u]>0)
    p2<-sum(newdata1[i,v]==0)
    print(i)
    if(p2==length(v) & p1>0){
      tmp<-c(i,p1)
      keep1<-rbind(keep1,tmp)
    }  
  }
  rlt$inforfreq<-keep1
  rlt$infornum<-nrow(keep1)
  
  Ofile2<-paste(dataset,"RF.informative.loci.freq.hist.in.cancer","jpeg",sep=".")
  jpeg(Ofile2)
  hist(keep1[,2]/length(u),main="",xlab="Hypermethylation Frequency in Cancer",lwd=1.5,cex.lab=1.5,cex.axis=1.5,col="red")
  dev.off()
  
  Ofile3<-paste(dataset,"informative.loci.importance.distribution.jpeg",sep=".")
  jpeg(Ofile3)
  plot(importance[order(importance[,2],decreasing=T)[1:num],2],main="",xlab="importance")
  dev.off()
  
  for(j in 1:num){
    a<-sum((importance[order(importance[,2],decreasing=T)[1:j],2])^2)
    b<-sum((importance[order(importance[,2],decreasing=T)[1:num],2])^2)
    if(a/b>0.8){
      bestnum<-j
      break
    } 
  }
  
  newdata2<-data[order(importance[,2],decreasing=T)[1:bestnum],]
  keep2<-c()
  for(i in 1:nrow(newdata2)){
    p1<-sum(newdata2[i,u]>0)
    p2<-sum(newdata2[i,v]==0)
    print(i)
    if(p2==length(v) & p1>0){
      tmp<-c(i,p1)
      keep2<-rbind(keep2,tmp)
    }  
  }
  rlt$top80freq<-keep2
  Ofile4<-paste(dataset,"RF.top80.informative.loci.freq.hist.in.cancer","pdf",sep=".")
  pdf(Ofile4)
  hist(keep2[,2]/length(u),main="",xlab="Hypermethylation Frequency in Cancer",lwd=1.5,cex.lab=1.5,cex.axis=1.5)
  dev.off()
  
  bestmtry1<-tuneRF(x=t(newdata2), y=as.factor(lab1), ntreeTry=5, stepFactor=2, improve=0.05,trace=TRUE, plot=TRUE, doBest=T)
  best.mtry1=bestmtry1$mtry
  rf.model4<-randomForest(t(newdata2),as.factor(lab1),mtry=best.mtry1,mtrees=500,importance=FALSE)
  rlt$acc2<-senspeacc(rf.model4)
  
  #importance2<-data.frame(ord=rownames(rf.model4$importance),importance=rf.model4$importance)
  colnames(newdata1)=colnames(file1)
  colnames(newdata2)=colnames(file1)
  Ofile5=paste(dataset,nrow(newdata1),"dataset","txt",sep=".")
  Ofile6=paste(dataset,nrow(newdata2),"dataset","txt",sep=".")
  write.table(newdata1,file=Ofile5,sep="\t",row.names=T,col.names=NA,quote=F)
  write.table(newdata2,file=Ofile6,sep="\t",row.names=T,col.names=NA,quote=F)
  
  return(rlt)
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
Rbedtoolsort<-function(functionstring="bedtools sort",bed1,opt.string=""){
  #create temp files
  a.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-i",a.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}


#######################################################################################################################
###   Analysis: 10/19/2015
#######################################################################################################################

# for desktop
setwd("C:/Users/shicheng/Dropbox/Project/methylation/monod")
sam<-read.table("sampleinfo.sort.txt",sep="\t",as.is=T)
new<-cbind(sapply(sam[,1],function(x) unlist(strsplit(x,":"))[2]),sam[,2])
new<-data.frame(new)
refgene<-read.table("hg19.bed")
write.table(new,file="sampleinfo.sort.new.txt",sep="\t",quote=F,col.names=NA,row.names=T)

# for server
setwd("/home/shg047/monod/mhl")
sam<-read.table("../sampleinfo.sort.txt",sep="\t",as.is=T)
new<-cbind(sapply(sam[,1],function(x) unlist(strsplit(x,":"))[2]),sam[,2])
new<-data.frame(new)
refgene<-read.table("../../annotation/hg19.bed")
write.table(new,file="../sampleinfo.sort.new.txt",sep="\t",quote=F,col.names=NA,row.names=T)

file1<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file2<-read.table("150209_BSPP_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file3<-read.table("140917_dRRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file4<-read.table("141216_SeqCap_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file5<-read.table("150209_SeqCap_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)

cor1<-match(colnames(file1),new[,1])
cor2<-match(colnames(file2),new[,1])
cor3<-match(colnames(file3),new[,1])
cor4<-match(colnames(file4),new[,1])
cor5<-match(colnames(file5),new[,1])

lab1<-new[cor1,2]
lab2<-new[cor2,2]
lab3<-new[cor3,2]
lab4<-new[cor4,2]
lab5<-new[cor5,2]

table(new[cor1,2])
table(new[cor2,2])
table(new[cor3,2])
table(new[cor4,2])
table(new[cor5,2])


rrbs<-pipeline(lab1=lab1,file1=file1,cor1=cor1,dataset="RRBS")
bspp<-pipeline(lab1=lab2,file1=file2,cor1=cor2,dataset="BSPP")
seqcap<-pipeline(lab1=lab5,file1=file5,cor1=cor5,dataset="SeqCap")

save(rrbs,file="rrbs.rlt.RData")
save(bspp,file="bspp.rlt.RData")
save(seqcap,file="seqcap.rlt.RData")

load(file="rrbs.rlt.RData")
load(file="bspp.rlt.RData")
load(file="seqcap.rlt.RData")


file1<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
samplename=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
cor1<-match(colnames(file1),new[,1])
lab1<-new[cor1,2]
groupname=lab1
colnames(matrix)=samplename
d <- dist(t(matrix)) # distance matrix
fit <- hclust(d, method="ward")         # distance matrix
pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.65,hang=-1,xlab="")
dev.off()


file2<-read.table("150209_BSPP_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
samplename=sapply(strsplit(colnames(file2),"[.]"),function(x) unlist(x)[1])
cor2<-match(colnames(file2),new[,1])
lab2<-new[cor2,2]
groupname=lab2
matrix=file2
colnames(matrix)=samplename
d <- dist(t(matrix),method="manhattan") # distance matrix
fit <- hclust(d, method="complete",xlab="")         # distance matrix
pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.85,hang=-1,xlab="BSPP Dataset Cluster Analysis")
dev.off()


file5<-read.table("150209_SeqCap_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
samplename=sapply(strsplit(colnames(file5),"[.]"),function(x) unlist(x)[1])
cor5<-match(colnames(file5),new[,1])
lab5<-new[cor5,2]
groupname=lab5
matrix=file5
colnames(matrix)=samplename
d <- dist(t(matrix),method="manhattan") # distance matrix
fit <- hclust(d, method="mcquitty")         # distance matrix
pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.85,hang=-1,xlab="SeqCap Dataset Cluster Analysis")
dev.off()




data<-t(file1)     
pheno=lab1
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

PCAPlot(t(file1),lab1,output="RRBS")
PCAPlot(t(file2),lab2,output="BSPP")
PCAPlot(t(file5),lab5,output="SeqCap")



head(file1)
barplot(data.matrix(file1))


A<-matrix(rnorm(1000,10,10),nrow=2,ncol=10)
boxplot(A)

plot(0,0,xlim=c(0,1),ylim=c(0,200),type="n",xlab="methylation level")
for(j in which(lab1=="Cancer")){
lines(density(file1[,j]))
j=j+1
}
for(j in which(lab1=="Normal")){
  lines(density(file1[,j]),col="blue")
}


cor<-cor(matrix,method="pearson")
sim<-1-cor
d <- dist(sim) # distance matrix
fit <- hclust(d, method="ward")         # distance matrix
pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.65)
dev.off()


rlt<-DataDescription(matrix=file1,colname=samplename,pheno=groupname)


library("gplots")
library("heatmap.plus")


list.files(pattern="*RData")

load("bspp.rlt.RData")
load("rrbs.rlt.RData")

u<-which(lab1=="Cancer")
v<-which(lab1=="Normal")
pvalue1<-apply(file1,1,function(x) t.test(x[u],x[v])$p.value)
length(which(pvalue1<0.05/length(pvalue1)))
heatmapdata<-file1[which(pvalue1<0.05/length(pvalue1)),]
head(heatmapdata)
heatmap.2(data.matrix(heatmapdata),col=redgreen(20))
heatmap.plus(data.matrix(heatmapdata))



u<-which(lab2=="Cancer")
v<-which(lab2=="Normal")
pvalue2<-apply(file2,1,function(x) t.test(x[u],x[v])$p.value)
length(which(pvalue2<0.05))

u1<-which(lab5=="Cancer")
v1<-which(lab5=="Normal")
pvalue5<-apply(file5,1,function(x) t.test(x[u],x[v])$p.value)
length(which(pvalue5<0.05/length(pvalue5)))


vcor<-cor(file1)
matrix=file1
DataDescription<-function(matrix,colname,pheno){
  ##matrix: column is individual, row is probe
  ##colname is Case1,case2,con1,con2 etc.
  ##pheno is ips,SCNT,ES and etc
  colnames(matrix)=colname
  ####### Figure1.dendrogram figure
  cor<-cor(matrix,method="pearson")
  disim<-1-cor
  d <- dist(disim) # distance matrix
  fit <- hclust(d, method="ward",scale=T)         # distance matrix
  save(fit,file="dendrogram.euclidean.ward.hclust.RData")
  pdf("Figure1.dendrogram.euclidean.ward.hclust.pdf")
  plot(fit,cex=0.5)
  dev.off()
  
  ####### Table1.correlation tables
  cor<-cor(matrix)
  write.table(cor,file="Table1.correlation.matrix.txt",sep="\t",quote=F,col.names=NA,row.names=T)
  
  ####### Figure2.PCA analysis
  data<-t(matrix)     
  pca <- prcomp(data,center=T,scale = F)  # Here, input file: row is individual and column is variable
  pdf("Figure2.PCA.loading.pdf")
  plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
  axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
  dev.off()
  var<-c()
  for(i in 1:length(pca$sdev)){
    var[i]<-sum((pca$sdev[1:i])^2)/sum((pca$sdev)^2)
  }
  pdf("Figure2.PCA.loading.total.increasing.pdf")
  plot(var,ylab="total variance",xlab="number of principle components",lwd=2)
  dev.off()
  save(pca,file="pca.result.RData")
  
  scores <- data.frame(pheno, pca$x[,1:3])
  col = as.numeric(as.factor(pheno))
  pdf("Figure2.PC12.pdf")
  plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.8,lwd=2)
  }
  dev.off()
  pdf("Figure2.PC13.pdf")
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  dev.off()
  pdf("Figure2.PC23.pdf")
  plot(x=scores$PC2,y=scores$PC3, xlim=c(min(scores$PC2),max(scores$PC2)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC2",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC2[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  dev.off()  
}





bedwithgap<-function(bed,gap){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
}

hybed1<-data.frame(cor2bed(rownames(rrbs$rawdata)[rrbs$hyperfreq[,1]]),rrbs$hyperfreq[,2])
hybed2<-data.frame(cor2bed(rownames(seqcap$rawdata)[seqcap$hyperfreq[,1]]),seqcap$hyperfreq[,2])
hybed3<-data.frame(cor2bed(rownames(bspp$rawdata)[bspp$hyperfreq[,1]]),bspp$hyperfreq[,2])


avg<-c()
for(gap in seq(1,300,by=1)){
gap=25
hybed1wgap<-bedwithgap(hybed1,gap)
hybed2wgap<-bedwithgap(hybed2,gap)
hybed3wgap<-bedwithgap(hybed3,gap)
head(rownames(rrbs$rawdata)[rrbs$hyperfreq[,1]])
head(hybed1)
head(hybed1wgap)
int1<-Rbedtools(bed1=hybed1wgap,bed2=hybed2wgap,opt.string="-wo")
nrow(int1)
int2<-Rbedtools(bed1=int1,bed2=hybed3wgap,opt.string="-wo")
nrow(int2)
bio<-paste(int2[,10],":",int2[,11]+gap,"-",int2[,12]-gap,sep="")
bsppDiag<-t(bspp$rawdata[match(bio,rownames(bspp$rawdata)),])
u<-which(lab2=="Cancer")
tmp<-c(gap,sum(rowSums(bsppDiag[u,])>0)/length(u))
avg<-rbind(avg,tmp)
print(tmp)
}

# then we can obtain the sensitivity trend along with the increasing of Gap
pdf("diagnosis.gap.pdf")
plot(avg[,1],avg[,2],type="l",xlab="Gap",ylab="Sensitivity",lwd=3,cex.axis=1.5,cex.lab=1.5,col="red")
dev.off()

#Obvously, Gap=25 is most perfect value (not very big and the sensivity is very large: Sen=93.75)   
gap=25
hybed1wgap<-bedwithgap(hybed1,gap)
hybed2wgap<-bedwithgap(hybed2,gap)
hybed3wgap<-bedwithgap(hybed3,gap)
head(rownames(rrbs$rawdata)[rrbs$hyperfreq[,1]])
head(hybed1)
head(hybed1wgap)
int1<-Rbedtools(bed1=hybed1wgap,bed2=hybed2wgap,opt.string="-wo")
nrow(int1)
int2<-Rbedtools(bed1=int1,bed2=hybed3wgap,opt.string="-wo")
nrow(int2)
bio<-paste(int2[,10],":",int2[,11]+gap,"-",int2[,12]-gap,sep="")
bsppDiag<-t(bspp$rawdata[match(bio,rownames(bspp$rawdata)),])
biomarker<-Rbedtools(bed1=int2,bed2=refgene,opt.string="-wo")
write.table(biomarker,file="biomarker.18.gene.txt",quote=F,col.names=F,row.names=F,sep="\t")
out<-unique(biomarker[,18])


# sensitivity increasing curve

input<-bsppDiag
mfvs<-function(data){
  rlt<-list()
  rlt$maxvariable<-which.max(apply(data,2,function(x) sum((x>0))))
  rlt$maxvalue<-max(apply(data,2,function(x) sum((x>0))))
  return(rlt)
}

value=1
sen<-c()
while(value>0){
  rlt<-mfvs(input)
  value=rlt$maxvalue
  sen<-c(sen,value)
  input<-input[-which(input[,rlt$maxvariable]>0),-rlt$maxvariable]
  print(c(value))
}


pdf("sen.93.7.pdf")
plot(cumsum(sen)/16,xlab="number",ylab="sensitivity",type="l",lwd=3,cex.lab=1.3)
dev.off()

x<-c(1:14)
sen<-c(4,2,2,1,1,1,1,1,1,1,0,0,0,0)
y<-cumsum(sen)/16
plot(x,y,type="o",lty=2,lwd=3,xlab="Number of biomarker",ylab="Sensitivity",col="blue",cex=1.5,cex.axis=1.5,cex.lab=1.5)
par<-par(mfrow=c(1,1))


# the range of the methylation frequency in the plasma is 6.25%-25%
freq<-apply(bsppDiag,2,function(x) sum(x[1:16]>0)/16)
range(freq)


range(int2[,4]/68)
range(int2[,8]/31)



bedfile2<-data.frame(cor2bed(rownames(file2)),file2)
int3<-Rbedtools(bed1=bedfile2,bed2=refgene,opt.string="-wao")
bedsort<-Rbedtoolsort(bed1=bedfile2)
data=bedsort

pdf("BSPP.distance.hist.pdf")
hist(uu,main="",xlab="Distance")
dev.off()


bedfile2<-data.frame(cor2bed(rownames(file2)),file2)
int3<-Rbedtools(bed1=bedfile2,bed2=bedfile2,opt.string="-wao")
bedsort<-Rbedtoolsort(bed1=bedfile2)
data=bedsort



out<-c()
rownames<-c()
for(i in 1:nrow(data)){
data<-Rbedtoolsort(bed1=data)
delta<-data[i+1,2]-data[i,3]
if(delta<gap){
  tmp<-data[i+1,4:ncol(data)]-data[i,4:ncol(data)]
  nametmp<-paste(data[i+1,1],":",data[i,2],"-",data[i+1,3],sep="")
  out<-rbind(out,tmp)
  rownames<-c(rownames,nametmp)
}else{
  out<-rbind(out,data[i,4:ncol(data)])
  rownames<-c(rownames,nametmp)
}
}


# for total loci
cor1<-rownames(file1)
cor2<-rownames(file2)
cor3<-rownames(file5)
bed1<-cor2bed(cor1)
bed2<-cor2bed(cor2)
bed3<-cor2bed(cor3)
int1<-Rbedtools(bed1=bed1,bed2=bed2)
int2<-Rbedtools(bed1=bed1,bed2=bed3)
int3<-Rbedtools(bed1=bed2,bed2=bed3)
int4<-Rbedtools(bed1=int1,bed2=bed3)
nrow(int1)
nrow(int2)
nrow(int3)
nrow(int4)
nrow(file1)
nrow(file2)
nrow(file3)

L1=as.numeric(as.character(bed1[,3]))-as.numeric(as.character(bed1[,2]))
L2=as.numeric(as.character(bed2[,3]))-as.numeric(as.character(bed2[,2]))
L3=as.numeric(as.character(bed3[,3]))-as.numeric(as.character(bed3[,2]))

quantile(L1)
quantile(L2)
quantile(L3)

# for informative loci
cor1<-rownames(file1)[rrbs$hyperfreq[,1]]
cor2<-rownames(file2)[bspp$hyperfreq[,1]]
cor3<-rownames(file5)[seqcap$hyperfreq[,1]]
bed1<-cor2bed(cor1)
bed2<-cor2bed(cor2)
bed3<-cor2bed(cor3)

int1<-Rbedtools(bed1=bed1,bed2=bed2,opt.string="-wo")
int2<-Rbedtools(bed1=bed1,bed2=bed3,opt.string="-wo")
int3<-Rbedtools(bed1=bed2,bed2=bed3,opt.string="-wo")
int4<-Rbedtools(bed1=int1,bed2=bed3,opt.string="-wo")
nrow(int1)
nrow(int2)
nrow(int3)
nrow(int4)

install.packages("VennDiagram")
library("VennDiagram")
# Make data
oneName <- function() paste(sample(LETTERS,5,replace=TRUE),collapse="")
geneNames <- replicate(1000, oneName())
GroupA <- sample(geneNames, 400, replace=FALSE)
GroupB <- sample(geneNames, 750, replace=FALSE)
GroupC <- sample(geneNames, 250, replace=FALSE)
GroupD <- sample(geneNames, 300, replace=FALSE)
v1 <- venn.diagram(list(A=GroupA, B=GroupB, C=GroupC), filename=NULL, fill=rainbow(3))
grid.newpage()
grid.draw(v1)


lab1=lab2
file1=file2
cor1=cor2
dataset="BSPP.2"











