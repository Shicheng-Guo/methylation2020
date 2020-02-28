setwd("/home/sguo/score")
i<-1
file=paste("chr",i,".funcPCA_score.txt",sep="")
data<-read.table(file,head=T,row.names=1,sep="\t",as.is=F)
for(i in 2:24){
  file=paste("chr",i,".funcPCA_score.txt",sep="")
  dat<-read.table(file,head=T,row.names=1,sep="\t",as.is=F)
  data<-cbind(data,dat)
}

write.table(data,"FunPCAScoreGenome.txt",sep="\t",col.names=NA, row.names=T,quote=F)

