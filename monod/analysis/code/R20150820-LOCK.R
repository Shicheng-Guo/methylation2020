
# download LDA from paper, liftover from hg18 to hg19 and save it to 512 server

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

setwd("C:\\Users\\shicheng\\Downloads")
setwd("/home/sguo/annotation")
list.files()
d2<-read.table("LOCK.hg19.bed",as.is=T,sep="\t")
d3<-read.table("LOCK.hg38.bed",as.is=T,sep="\t")
head(d2)
head(d3)
bed2<-data.frame(cor2bed(d2[,1]),d2[,1])
bed3<-data.frame(cor2bed(d3[,1]),d3[,1])
head(bed2)
head(bed3)

write.table(bed2,file="LOCK.hg19.bed4",col.names=F,row.names=F,quote=F,sep="\t")
write.table(bed3,file="LOCK.hg38.bed4",col.names=F,row.names=F,quote=F,sep="\t")




