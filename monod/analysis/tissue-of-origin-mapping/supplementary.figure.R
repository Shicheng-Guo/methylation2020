
data<-read.table("C:\\Users\\shicheng\\Downloads\\RASSF1.chr3.rsq",head=T,sep="\t",row.names = 1)
data[1:5,1:5]
library("grDevices")
col=colorRampPalette(c("white", "red"))(20)
M<-data.matrix(data)
M[lower.tri(M)] <- NA
M[is.na(M)]<-0
pdf("RASSF1.pdf")
image(M,col = col,frame=F,xaxt="n",yaxt="n")
dev.off()


data<-read.table("C:\\Users\\shicheng\\Downloads\\SHOX2.chr3.rsq",head=T,sep="\t",row.names = 1)
data[1:5,1:5]
library("grDevices")
col=colorRampPalette(c("white", "red"))(20)
M<-data.matrix(data)
M[lower.tri(M)] <- NA
M[is.na(M)]<-0
pdf("SHOX2.pdf")
image(M,col = col,frame=F,xaxt="n",yaxt="n")
dev.off()

data<-read.table("RASSF1.61-WGBS.coverage",head=F,sep="\t")
coverage<-data[,6]

r2<-read.table("/home/shg047/work/monod/Fig1B/hapinfo/RASSF1.chr3.rsq",sep="\t",row.names=1,head=T,check.names = F)
coverage<-read.table("/home/shg047/work/monod/Fig1B/bam/RASSF1.61-WGBS.coverage",sep="\t")
data<-data.frame(pos=coverage[,5]+coverage[,2],cov=coverage[,6])
input<-data[match(rownames(r2),data[,1]),2]
pdf("RASSF1.pdf")
barplot(input,col="blue",border=F)
dev.off()


setwd("C:\\Users\\shicheng\\Downloads\\")

data<-read.table("C:\\Users\\shicheng\\Downloads\\scnt.wnt5a.chr14.rsq",head=T,sep="\t",row.names = 1)
data[1:5,1:5]
dim(data)
library("grDevices")
col=colorRampPalette(c("white", "red"))(20)
M<-data.matrix(data)
M[lower.tri(M)] <- NA
M[is.na(M)]<-0
pdf("scnt.wnt5a.chr14.rsq.pdf")
image(M,col = col,frame=F,xaxt="n",yaxt="n")
dev.off()

data<-read.table("C:\\Users\\shicheng\\Downloads\\ipsc.wnt5a.chr14.rsq",head=T,sep="\t",row.names = 1)
data[1:5,1:5]
library("grDevices")
col=colorRampPalette(c("white", "red"))(20)
M<-data.matrix(data)
M[lower.tri(M)] <- NA
M[is.na(M)]<-0
pdf("ipsc.wnt5a.chr14.rsq.pdf")
image(M,col = col,frame=F,xaxt="n",yaxt="n")
dev.off()


coverage<-read.table("C:\\Users\\shicheng\\Downloads\\ipsc.wnt5a.cov",sep="\t")
cov<-data.frame(pos=coverage[,5]+coverage[,2]-1,cov=coverage[,6])
input<-cov[match(rownames(data),cov[,1]),2]
length(input)
pdf("ipsc.wnt5a.cov.pdf")
barplot(input,col="blue",border=F)
dev.off()

coverage<-read.table("C:\\Users\\shicheng\\Downloads\\scnt.wnt5a.cov",sep="\t")
cov<-data.frame(pos=coverage[,5]+coverage[,2]-1,cov=coverage[,6])
input<-cov[match(rownames(data),cov[,1]),2]
length(input)
pdf("scnt.wnt5a.cov.pdf")
barplot(input,col="blue",border=F)
dev.off()


