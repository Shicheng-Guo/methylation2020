# this code was used to check whether the distribnution of MHL would change for WGBS dataset with the regions of RRBS

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
  return(cor)
}

setwd("/home/shg047/monod/phase2")
file1<-read.table("WGBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file2<-read.table("RRBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
bed1<-cor2bed(rownames(file1))
bed2<-cor2bed(rownames(file2))
bed<-Rbedtools(functionstring="intersectBed",bed1,bed2,opt.string="-wa -u")
cor1<-gsub("[ ]","",bed2cor(bed))

file1<-file1[match(cor1,rownames(file1)),]
colnames(file1)
colnames(file1)<-gsub("_","-",colnames(file1))
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
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
matrix<-matrix[,-c(11,12)]
d <- dist(t(matrix)) # distance matrix
fit <- hclust(d, method="complete")         # distance matrix

setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\monod")
load("WGBS.RRBS.RData")
par(mar=c(2,7,1,1))
boxplot(matrix[,1:ncol(matrix)],outline=F,horizontal=T,notch=T,las=1,cex.axis=0.65)


