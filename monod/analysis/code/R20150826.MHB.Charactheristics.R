#####################################################################
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com \
###   updata time: 5/29/2015
###   Figure 1. Methylation Dataset Collection 
###   Figure 2. PCA Analysis
###   Fiugre 3. Differential Analysis
###   Figure 4. Pathway Analysis
###   Figure 5. heatmap
#####################################################################

file="WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed"
cgi="/home/shg047/bioin/annotation/hg19/hg19.cpgIslandExt.txt"
cgshore="/home/shg047/bioin/annotation/hg19/hg19.cpgshoreExt.bed"
cgshelf="/home/shg047/bioin/annotation/hg19/hg19.cpgshelfExt.bed"
refgene="/home/shg047/bioin/annotation/hg19/hg19_refGene.bed"
repeatmask="/home/shg047/bioin/annotation/hg19/hg19.repeatmasker.bed"
tss="/home/shg047/bioin/annotation/hg19/hg19.refGeneTSS.bed"

wc -l $file
wc -l $cgi
wc -l $cgshore
wc -l $cgshelf
wc -l $refgene


bedtools intersect -u -a $file -b $cgi | wc -l 
bedtools intersect -u -a $file -b $cgshore | wc -l 
bedtools intersect -u -a $file -b $cgshelf | wc -l 
bedtools intersect -u -a $file -b $refgene | wc -l 
bedtools intersect -wb  -a $file -b $repeatmask | wc -l 
bedtools closest -D -a $file -b $tss | head

bedtools closest -D b -a WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed -b /home/shg047/bioin/annotation/hg19/hg19.refGeneTSS.bed > WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed.Tss.Dis

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\phase3")
mhb<-read.table("WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed")
head(mhb)
len<-mhb[,3]-mhb[,2]
hist(len[len<350],breaks=100,main="",xlab="",ylab="",col="blue")



R<-read.table("WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed.Tss.Dis")
hist(R[,ncol(R)],xlim=c(-100000,100000),main="",breaks=750,xlab="",cex.lab=1.15,cex.axis=1.15,lwd=1.15,col="blue",col.axis="blue",col.lab="blue")

wc -l $file
wc -l $refgene
wc -l $repeatmask

bedtools intersect -u -a $file -b $refgene | wc -l
bedtools intersect -u -a $file -b $repeatmask | wc -l
bedtools intersect -u -a $refgene -b $repeatmask | wc -l 
bedtools intersect -u -a $file -b $refgene | bedtools intersect -u -a - -b $repeatmask | wc -l

library("VennDiagram")
# Venn for PBMC
area1=49542
area2=789065
area3=5232241
n12=34095
n13=20714
n23=367983
n123=13017
jpeg("Figuare 1c.characteristic of methylation blocks in human genome.jpg")
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, cex=2,cat.cex=1.5,col.lab="white",cat.pos=c(350,0,180),category =c("Methylation Block","RefGene","Repeat"),col=c(2:4),fill=2:4,lwd=2,ind = T, list.order = 1:3)
dev.off()


zcat file.gz | less
less file.gz; 
zcat $file.gz | cut -f 1-10 | less



d<-read.table("WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed.maskrepeat")
table(d[,10])

q<-c(20714,5845,4269,2313,34095,15447)/49542
q



setwd("/home/sguo/monod/phase3")
wgbsmb<-read.table("WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",sep="\t")
load("/home/sguo/monod/rrbs/RRBS.MethylationBlock.RData")
rrbs<-data.frame(rrbsmb$hdr,rrbsmb$hdrc)
rrbsmb<-subset(rrbs,ratio>0.6)
rrbsmbc<-data.frame(cor2bed(rownames(rrbsmb)),rownames(rrbsmb))
write.table(rrbsmbc,file="rrbs.encode.mb.0.6.bed",sep="\t",col.names=F,row.names=F,quote=F)

out<-c()
for(i in seq(0.2,0.999,by=0.05)){
  sub<-subset(rrbs,ratio>i)
  bed<-cor2bed(rownames(sub))
  tmp<-Rbedtools("intersectBed",bed,wgbsmb,opt.string="-wa")
  out<-rbind(out,c(i,nrow(tmp),nrow(tmp)/nrow(sub)))
}

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\RRBS")
load("rrbs.wgbs.mb.overlap.RData")
plot(x=out[,1],y=out[,3])
qplot(x=out[,1],y=out[,3])


#!/bin/bash
for i in {1..1000}
do
OBS=`bedtools shuffle -i rrbs.encode.mb.0.6.bed -g /home/sguo/annotation/hg19.chrom.sizes | bedtools intersect -a WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed -b - | wc -l`
echo $OBS 
done

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\phase3")



#####################################################################
###################### Functions and librarys  ######################
#####################################################################

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
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


