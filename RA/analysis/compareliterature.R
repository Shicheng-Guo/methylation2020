# compare our DMS with others DMS
setwd("/home/shg047/meth450/dyh/analysis")
dms1<-read.table("sdms.txt",sep="\t",head=T,as.is=T)
dms2<-read.table("PMC3549371.txt",sep="\t",as.is=T)
dms3<-read.table("PMID24005183.txt",sep="\t",head=T,as.is=T)
dms4<-read.table("eGwas.txt",sep="\t",head=F,as.is=T)



df1<-dms1[na.omit(match(dms2$V5,dms1$V5)),]
df2<-dms1[na.omit(match(dms3$Gene.symbol,dms1$V5)),]
df3<-dms3[na.omit(match(dms2$V5,dms3$Gene.symbol)),]
df4<-dms1[na.omit(match(dms4[,1],dms1[,2])),]

dms2[match("GALNT9",dms2$V5),]