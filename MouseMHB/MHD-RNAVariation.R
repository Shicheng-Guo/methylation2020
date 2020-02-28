# This software is Copyright © 2017 The Regents of the University of California. All Rights Reserved.
#  
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
#  
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# kzhang@ucsd.edu

##===============================================================================================================================================================
##   Try to identify the relationship between methylation linkage disequlibrium and gene expression variation in iPS and SCNT cells. 
##
##================================================================================================================================================================
library("tidyr")
library("reshape2")
library("ggplot2")
library("vioplot")
##===============================================================================================================================================================
##   Prepare Database and shared annoation 
##   house keeping gene (hpk)
##================================================================================================================================================================
hkg<-read.table("/media/Home_Raid1/shg047/work/Alice/mouse/housekeeping.txt",sep="\t")
exp_GSE75804<-read.table("/media/NAS3_volume2/shg047/Alice/mouse/entropy/GSE75804_combined.counts.txt",sep="\t",row.names=1,head=T)
plogenes<-read.table("/home/shg047/work/Alice/mouse/plogene.txt") # first column gene name
##===============================================================================================================================================================
##   Load functions
##   Clust analysis, PCA, Prediction
##================================================================================================================================================================

# setwd("/media/Home_Raid1/zhl002/WGBS_mouse")
cp /media/Home_Raid1/zhl002/WGBS_mouse/allsample.txt /media/Home_Raid1/shg047/NAS3/Alice/mouse/entropy
cp /media/Home_Raid1/zhl002/WGBS_mouse/iPSC_SCNT.R2Pvalue.high.txt /media/Home_Raid1/shg047/NAS3/Alice/mouse/entropy

# reference
hkg<-read.table("/media/Home_Raid1/shg047/work/Alice/mouse/housekeeping.txt",sep="\t")
exp_GSE75804<-read.table("/media/NAS3_volume2/shg047/Alice/mouse/entropy/GSE75804_combined.counts.txt",sep="\t",row.names=1,head=T)

# plo plogenes 

hkfreq(exp_GSE75804,output="hkgdatatable-GSE75804.txt")

hkfreq<-function(file,output="hkgdatatable-GSE75804.txt"){
hkgdata<-file[na.omit(match(as.character(hkg$V4),rownames(file))),]
hkgdatatable<-apply(hkgdata,1,function(x) sum(x>0)/ncol(hkgdata))
write.table(hkgdatatable,file=output,sep="\t",col.names=NA,row.names=T,quote=F)
return(hkgdatatable)
}

setwd("/media/Home_Raid1/shg047/NAS3/Alice/mouse/entropy")

file1<-read.table("/media/Home_Raid1/zhl002/WGBS_mouse/allsample.txt")
file2<-read.table("/media/Home_Raid1/zhl002/WGBS_mouse/iPSC_SCNT.R2Pvalue.high.txt")
refGene<-read.table("/media/Home_Raid1/shg047/work/db/mm9/mm9.refGene.bed")
head(refGene)
# bedtools intersect -a bed1.bed -b /media/Home_Raid1/shg047/work/db/mm9/mm9.refGene.bed -wao > xm.txt
bed1<-cor2bed(as.character(file2[,1]))
bed2<-refGene
# write.table(bed1,file="bed1.bed",sep="\t",quote=F,row.names=F,col.names=F)
bed<-bed2gene(bed1,bed2)
bed1<-subset(bed,V10=="Promoter")
bed2<-subset(bed,V10=="Enhancer")
bed3<-subset(bed,V10=="Exon")
bed4<-subset(bed,V10=="Intron")
bed5<-subset(bed,V10=="UTR5")
bed6<-subset(bed,V10=="UTR3")

####################### while regions(enhancer,promoter,exon,intron) ####################### 
input<-file1
# calculate sd or entropy for groups and then plot figures
x1<-grep("ips",colnames(input))
x2<-grep("NT",colnames(input))

iPS<-apply(input,1,function(x){ x = discretize(x, numBins=10);shannon(x[x1])})
SCNT<-apply(input,1,function(x){ x = discretize(x, numBins=10);shannon(x[x2])})

p = discretize(input, numBins=10)

y=data.frame(scale(cbind(iPS,SCNT)))
library("preprocessCore")
yqn<-normalize.quantiles(data.matrix(y))
rownames(yqn)=rownames(y)
yqnsubset<-yqn[unique(na.omit(match(subset(bed,V10=="Promoter"|V10=="Enhancer")$V9,rownames(yqn)))),]
colnames(yqnsubset)<-c("iPS","SCNT")
wilcox.test(yqnsubset[,1],yqnsubset[,2],paired = T)

input.long<-melt(y)
head(input.long)
input.long <- within(input.long,variable=factor(variable,levels=c("iPS","SCNT")))
pdf("entroy.full-region.boxplot.whole.pdf")
ggplot(aes(y = value, x = variable), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw()+ 
  geom_point(position = position_jitter(width = 0.2))+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
  )
dev.off()
pdf("entroy.full-region.violin.whole.pdf")
vioplot(iPS, SCNT, names=c("iPS", "SCNT"), col="gold")
dev.off()

# calculate sd or entropy for groups and then plot figures
y2<-unique(na.omit(match(bed1$V9,rownames(y))))
y3<-unique(na.omit(match(bed2$V9,rownames(y))))
y4<-unique(na.omit(match(bed3$V9,rownames(y))))
y5<-unique(na.omit(match(bed4$V9,rownames(y))))
y6<-unique(na.omit(match(bed5$V9,rownames(y))))

pyp(y2)
pyp(y3)
pyp(y4)
pyp(y5)
pyp(y6)

pyp<-function(yk){
yp<-y[yk,]
p<-t.test(yp[,1],yp[,2],paired = T)$p.value
return(p)
}

wilcox.test(y2[,1],y2[,2],paired = T)
input.long<-melt(y,id.vars = c("iPS", "SCNT"))
input.long <- within(input.long,Var2 <- factor(Var2,levels=c("iPS","SCNT")))
pdf("entroy.ips-high-R2.boxplot.whole.pdf")
ggplot(aes(y = value, x = Var2), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw()+ 
  geom_point(position = position_jitter(width = 0.2))+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
  )
dev.off()
pdf("entroy.ips-high-R2.violin.whole.pdf")
vioplot(iPS, SCNT, names=c("iPS", "SCNT"), col="gold")
dev.off()


####################### promoter and enhancer ####################### 
# promoter and enhancer
input<-file1[unique(na.omit(match(subset(bed,V10=="Promoter"|V10=="Enhancer")$V9,rownames(file1)))),]
# calculate sd or entropy for groups and then plot figures
library("entropy")
x1<-grep("ips",colnames(input))
x2<-grep("NT",colnames(input))
shannon(as.numeric(data.matrix(input[,x1])))
shannon(as.numeric(data.matrix(input[,x2])))
iPS<-apply(input[,x1],1,function(x) shannon(x))
SCNT<-apply(input[,x2],1,function(x) shannon(x))
y=cbind(iPS,SCNT)
wilcox.test(y[,1],y[,2],paired = T)
input.long<-melt(y,id.vars = c("iPS", "SCNT"))
input.long <- within(input.long,Var2 <- factor(Var2,levels=c("iPS","SCNT")))

pdf("entroy.ips-high-R2.boxplot.pdf")
ggplot(aes(y = value, x = Var2), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw()+ 
  geom_point(position = position_jitter(width = 0.2))+
  theme(plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank()
)
dev.off()

# promoter and enhancer
input<-file1[unique(na.omit(match(subset(bed,V10=="Intron"|V10=="Exon")$V9,rownames(file1)))),]
# calculate sd or entropy for groups and then plot figures
library("entropy")
x1<-grep("ips",colnames(input))
x2<-grep("NT",colnames(input))
shannon(as.numeric(data.matrix(input[,x1])))
shannon(as.numeric(data.matrix(input[,x2])))
iPS<-apply(input[,x1],1,function(x) shannon(x))
SCNT<-apply(input[,x2],1,function(x) shannon(x))
y=cbind(iPS,SCNT)
wilcox.test(y1,y2,paired = T)
input.long<-melt(y,id.vars = c("iPS", "SCNT"))
input.long <- within(input.long,Var2 <- factor(Var2,levels=c("iPS","SCNT")))

pdf("entroy.ips-high-R2.boxplot-exon-intron.pdf")
ggplot(aes(y = value, x = Var2), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw()+ 
  geom_point(position = position_jitter(width = 0.2))+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
  )
dev.off()



###### Function loaded before Analysis ##########
shannon <- function(p){
  library("entropy")
  if (min(p) < 0 || sum(p) <= 0) return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

bed2gene<-function(bed1,bed2){
  bed<-Rbedtools(functionstring="intersectBed",bed1,bed2,opt.string="-wao")
  return(bed)
}

bedwithgap<-function(bed,gap=2000){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
}

cor2bed<-function(cor){
  cor<-as.character(cor)
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

bed2cor<-function(bed){
  bed<-data.frame(bed)
  a<-unlist(apply(bed,1,function(x) paste(x[1],":",x[2],"-",x[3],sep="")))
  cor<-as.character(a)
  return(cor)
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





