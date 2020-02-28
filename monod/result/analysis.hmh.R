par(mfrow=c(2,2))


d<-read.table("CRC-MHM4.txt",sep="\t",head=T,check.names = F,as.is=T)
bed<-d[,1]
d<-data.matrix(d[,2:ncol(d)])
d[1:5,1:5]
id1<-grep("NC-|STL",colnames(d))
id2<-grep("6-P",colnames(d))

iid<-which(apply(d,1,function(x) sum(x[id1]))==0)

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

write.bed<-function(bed,file,extend=0){
  bed[,2]<-as.numeric(as.character(bed[,2]))-extend
  bed[,3]<-as.numeric(as.character(bed[,3]))+extend
  if(ncol(bed)==3){
  bed[,4]<-paste(bed[,1],":",bed[,2],"-",bed[,3],sep="")  
  }
  if(ncol(bed)>=4){
    write.table(bed,file=file,sep="\t",col.names=F,row.names=F,quote=F)
  }
}

write.bed(cor2bed(bed[iid]),"HMH-6.bed")


dmax<-d[iid,id2]
dmax[dmax>1]=1
hist(apply(dmax,1,sum),xlim=c(0,20),main="CRC",xlab="Sample Number")
iid[which.max(apply(dmax,1,sum))]


d<-read.table("C:\\Users\\shicheng\\Downloads\\HMH-7.txt",sep="\t",head=T,check.names = F,as.is=T)
bed<-d[,1]
d<-data.matrix(d[,2:ncol(d)])
id1<-grep("NC-|STL",colnames(d))
id2<-grep("7-P",colnames(d))

iid<-which(apply(d,1,function(x) sum(x[id1]))==0)
write.bed(cor2bed(bed[iid]),"HMH-7.bed")



dd1<-read.table("CRC.High.DMHL.txt",head=F)
write.bed(cor2bed(as.character(dd1[,1])),"CRC.DMHL.bed",5000)
dd2<-read.table("LC.High.DMHL.txt",head=F)
write.bed(cor2bed(as.character(dd2[,1])),"LC.DMHL.bed",5000)

bedtools intersect -wa -a CRC.DMHL.bed -b CRC-MHM4.txt | wc -l 
bedtools intersect -wa -a LC.DMHL.bed -b LC-MHM4.txt  | wc -l 

bedtools intersect -wa -a CRC.DMHL.bed -b CRC-MHM4.txt > CRC.Share.bed
bedtools intersect -wa -a LC.DMHL.bed -b LC-MHM4.txt  > LC.Share.bed

bedtools intersect -wa -b CRC.Share.bed -a ~/db/hg19/hg19_refGene.bed
bedtools intersect -wa -b LC.Share.bed -a ~/db/hg19/hg19_refGene.bed


x<-c(0.1,0.2,0.3,0.4)
y<--2*x
rlt<-cor(x,y)
rlt
cov(x,y)
? cor

x<-1:100
barplot(x,col="green",border = F)
? barplot

d<-read.table("CRC-MHM4.txt",sep="\t",head=T,check.names = F,as.is=T)






dmax<-d[iid,id2]
dmax[dmax>1]=1
hist(apply(dmax,1,sum),xlim=c(0,20),main="CRC",xlab="Sample Number")
iid[which.max(apply(dmax,1,sum))]


dmax<-d[which(apply(d,1,function(x) sum(x[id1]))==0),id2]
dmax[dmax>1]=1
hist(apply(dmax,1,sum),xlim=c(0,20),main="LC",xlab="Sample Number")
max(apply(dmax,1,sum))



d<-read.table("C:\\Users\\shicheng\\Downloads\\data.input.txt")
IQR(d[,1])
