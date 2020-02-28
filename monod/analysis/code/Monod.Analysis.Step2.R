setwd("/home/shg047/monod/phase2")

dat1<-read.table("RRBS_methHap_load_matrix_July2015.txt",head=T,row.names=1,sep="\t",as.is=T)
dat2<-read.table("WGBS_methHap_load_matrix_July2015.txt",head=T,row.names=1,sep="\t",as.is=T)
dat3<-read.table("WGBS_SeqCap_methHap_load_matrix_July2015.txt",head=T,row.names=1,sep="\t",as.is=T)

setwd("/home/shg047/monod/phase2")
sam<-read.table("/home/shg047/monod/sampleinfo.sort.txt",sep="\t",as.is=T)
sam[,1]<-gsub("6P","6-P",sam[,1])
sam[,1]<-gsub("7P","7-P",sam[,1])
sam[,1]<-gsub("6T","6-T",sam[,1])
sam[,1]<-gsub("7T","7-T",sam[,1])
sam[,1]<-gsub("NCP","NC-P",sam[,1])
sam[,1]<-gsub("NCT","NC-T",sam[,1])
sam[,1]<-gsub("PCP","PC-P",sam[,1])
sam[,1]<-gsub("PCT","PC-T",sam[,1])
sam[,1]<-gsub("_1","",sam[,1])
sam[,1]<-gsub("_2","",sam[,1])
new1<-as.character(sapply(sam[,1],function(x) unlist(strsplit(x,":"))[2]))
new2<-as.character(sapply(new1,function(x) unlist(strsplit(x,"[.]"))[1]))
new<-data.frame(cbind(sam[,1],new1,new2,sam[,2]))
write.table(new,file="saminfo.txt",col.names=F,row.names=F,quote=F,sep="\t")

file1<-read.table("RRBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
samplename2<-gsub("NC-P-12","NC-12",samplename2)

cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
matrix=file1
colnames(matrix)=samplename

d <- dist(t(matrix)) # distance matrix
fit <- hclust(d, method="ward")         # distance matrix
pdf("Figure1.dendrogram.pearson.ward.hclust.pdf")
plot(fit,cex=0.45,hang=-1,xlab="")
dev.off()


file2<-read.table("WGBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
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


file5<-read.table("WGBS_SeqCap_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
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
