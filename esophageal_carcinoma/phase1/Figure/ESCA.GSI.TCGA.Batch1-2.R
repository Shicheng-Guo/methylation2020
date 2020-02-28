
# Batch 1: 82 region. However, company found large number of them can not designed out the primer. 

setwd("/home/shg047/oasis/TCGA/Meth/Data")
data<-read.table("./split/MergePancancerMethMatrix.txt",sep="\t",as.is=T)
candidate<-subset(data,data[,2]=="ESCA_1" & data[,3]>0.45)
dim(candidate)
anno<-read.table("~/oasis/db/GPL13534.sort.bed",sep="\t",as.is=T)
rlt<-anno[match(candidate[,1],anno[,4]),]
table(rlt[,5])
library("ggplot2")
args<-c("header.txt","PancancerMethSaminfo_March2016.txt")
colname<-read.table(args[1],sep="\t",row.names=1,colClasses=c("character",rep("character",6440)),nrow=1)
saminfo<-read.table(args[2],sep="\t")
group1<-paste(saminfo[match(colname,saminfo[,1]),3],"_",saminfo[match(colname,saminfo[,1]),4],sep="")
group2<-saminfo[match(colname,saminfo[,1]),3]
group3<-saminfo[match(colname,saminfo[,1]),4]
for(i in 1:nrow(candidate)){
cpg<-candidate[i,1]
cmd<-paste("grep ",cpg," PancancerMethMatrix_March2016.txt > ",cpg,".txt",sep="")
system(cmd)
data<-read.table(paste(cpg,".txt",sep=""),sep="\t",colClasses=c("character",rep("numeric",6440)),nrow=2500,row.names=1)
newdata=data.frame(beta=as.numeric(data),SamType=group3,cancerType=group2,CancerType=group1)
newdata<-na.omit(newdata)
ggplot(newdata, aes(factor(CancerType),beta)) + geom_boxplot(aes(fill = factor(SamType)))+ coord_flip()
outputFile=paste(cpg,".pdf",sep="")
ggsave(outputFile)
}


# Batch 2: 82 region. However, company found large number of them can not designed out the primer. 
setwd("/home/shg047/oasis/TCGA/Meth/Data")
data<-read.table("./split/MergePancancerMethMatrix.txt",sep="\t",as.is=T)
candidate<-subset(data,data[,2]=="ESCA_1" & data[,3]<=0.45 & data[,3]>=0.38)
dim(candidate)
anno<-read.table("~/oasis/db/GPL13534.sort.bed",sep="\t",as.is=T)
rlt<-cbind(candidate,anno[match(candidate[,1],anno[,4]),])
table(rlt[,5])
length(table(rlt[,5]))
write.table(rlt, file="Batch-3.ESCA.high.GSI.CpGs.txt",col.names=NA,row.names=T,sep="\t",quote=F)

library("ggplot2")
args<-c("header.txt","PancancerMethSaminfo_March2016.txt")
colname<-read.table(args[1],sep="\t",row.names=1,colClasses=c("character",rep("character",6440)),nrow=1)
saminfo<-read.table(args[2],sep="\t")
group1<-paste(saminfo[match(colname,saminfo[,1]),3],"_",saminfo[match(colname,saminfo[,1]),4],sep="")
group2<-saminfo[match(colname,saminfo[,1]),3]
group3<-saminfo[match(colname,saminfo[,1]),4]
for(i in 1:nrow(candidate)){
  cpg<-candidate[i,1]
  cmd<-paste("grep ",cpg," PancancerMethMatrix_March2016.txt > ",cpg,".txt",sep="")
  system(cmd)
  data<-read.table(paste(cpg,".txt",sep=""),sep="\t",colClasses=c("character",rep("numeric",6440)),nrow=2500,row.names=1)
  newdata=data.frame(beta=as.numeric(data),SamType=group3,cancerType=group2,CancerType=group1)
  newdata<-na.omit(newdata)
  ggplot(newdata, aes(factor(CancerType),beta)) + geom_boxplot(aes(fill = factor(SamType)))+ coord_flip()
  outputFile=paste(cpg,".pdf",sep="")
  ggsave(outputFile)
  print(i)
}

# Step 3. remove PBMC hyper-CpGs.
pbmc<-read.table("/home/shg047/oasis/TCGA/Meth/Normal.PBMC.GEO.HM450K.Beta.txt",head=T,row.names=1,sep="\t",as.is=T)
rlt2<-pbmc[match(rlt[,4],rownames(pbmc)),]
write.table(rlt2,file="batch3.CpGs.in.PBMC.txt",col.names=NA,row.names=T,sep="\t",quote=F)



