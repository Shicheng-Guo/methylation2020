

setwd("/home/shg047/oasis/TCGA/Meth/Data")
candidate<-read.table("weilin.txt",sep="\t",as.is=T)
library("ggplot2")
args<-c("header.txt","PancancerMethSaminfo_March2016.txt")
colname<-read.table(args[1],sep="\t",row.names=1,colClasses=c("character",rep("character",6440)),nrow=1)
saminfo<-read.table(args[2],sep="\t")
group1<-paste(saminfo[match(colname,saminfo[,1]),3],"_",saminfo[match(colname,saminfo[,1]),4],sep="")
group2<-saminfo[match(colname,saminfo[,1]),3]
group3<-saminfo[match(colname,saminfo[,1]),4]
for(i in 1:nrow(candidate)){
  cpg<-candidate[i,4]
  cmd<-paste("grep ",cpg," PancancerMethMatrix_March2016.txt > ",cpg,".txt",sep="")
  system(cmd)
  data<-read.table(paste(cpg,".txt",sep=""),sep="\t",colClasses=c("character",rep("numeric",6440)),nrow=2500,row.names=1)
  newdata=data.frame(beta=as.numeric(data),SamType=group3,cancerType=group2,CancerType=group1)
  newdata<-na.omit(newdata)
  ggplot(newdata, aes(factor(CancerType),beta)) + geom_boxplot(aes(fill = factor(SamType)))+ coord_flip()
  outputFile=paste(cpg,".pdf",sep="")
  ggsave(outputFile)
}


