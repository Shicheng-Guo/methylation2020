# source("http://bioconductor.org/biocLite.R")
# biocLite("quantsmooth")

library("quantsmooth")
setwd("/home/shg047/monod/methyblock")
file<-"WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed";
data<-read.table(file,head=F,sep="\t",as.is=T)
colnames(data)=c("CHR","Pos","End","Length")
head(data)



data<-read.table("C:\\Users\\shicheng\\Downloads\\chrom.plot.txt",head=T,row.names=1,sep="\t",as.is=T)
chrompos<-data
data[,1]<-gsub("chr","",data[,1])
pdf("chrosomeplot.pdf")
par(mar=c(0.1,1,1,1))
CHR<-data[,1]  # Chromosomes
MapInfo<-data[,2]# position on chromosome
chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE, organism="hsa",sexChromosomes = TRUE,units="bases",cols="green")
points(chrompos[,2],chrompos[,1]+0.2,pch="|",col="red",cex=1.1)
dev.off()


