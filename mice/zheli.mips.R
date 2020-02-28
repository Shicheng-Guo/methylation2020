# Genome-wide DNA methylation of miPS and 
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


setwd("/media/LTS_33T/ZL_LTS33/mouse_WGBS/bedfiles")

pvalue<-read.table("output_DMS.tsv")
perl -lane 'print "@F[0]\t@F[1]\t@F[0]\t@F[0]\t@F[0]\t@F[0]\t"'
sig<-read.table("sig.txt",head=T,sep="\t",as.is=T)
pos2=sig$pos+1
refgene<-read.table("/home/shg047/annotation/hg19.bed")
newsig<-data.frame(sig$chr,sig$pos,pos2,sig$p)
Rbedtools(bed1=newsig,bed2=refgene,opt="-wao")

pdf("delta.dmr.pdf")
hist(pvalue$diff,col="red",breaks=50)
dev.off()


# pvalue calculation

data<-read.table("pvalue.shore.txt")
number<-34
p=1-sum(number>data[,1])/nrow(data)

  
for i in `ls *.BED.txt`
do 
bedtools intersect -a DSS.DMR.Shore.34.bed -b $i -wb > $i.shore 
done
  

file<-list.files(pattern="*.txt.shore")
cor<-c()
for(i in 1:length(file)){
  f=read.table(file[i],sep="\t",as.is=T)
  u<-paste(f[,1],":",f[,2],"-",f[,3],sep="")
  v<-data.frame(u,f)
  cor<-c(cor,u)
}

matrix<-c()
for(i in 1:length(file)){
  f=read.table(file[i],sep="\t",as.is=T)
  u<-paste(f[,1],":",f[,2],"-",f[,3],sep="")
  tmp<-f[match(cor,u),7]
  matrix<-cbind(matrix,tmp)
}
colnames(matrix)<-file
rownames(matrix)<-cor

write.table(matrix,file="matrix.detail.mice.alicey.txt",col.names=NA,row.names=T,sep="\t",quote=F)

# collect the data of 
for i in `ls GSM*.tsv`
do
awk '$6>=5 && $7==1' $i > $i.trim &
done


# get the sample ID
ls *trim | awk -F_ '{print $1}' | sort -u 


# loop the sample id and print them
for i in GSM1385973 GSM1385974 GSM1385975 GSM1385976 GSM1385977 GSM1385978 GSM1385979 GSM1385980 GSM1385981 GSM1385982 GSM1385983
do 
for j in `ls $i*.trim`
do
cat $j | sed 1d >> $i.dat
done
done

(2/21614)


for i in GSM1385973 GSM1385974 GSM1385975 GSM1385976 GSM1385977 GSM1385978 GSM1385979 GSM1385980 GSM1385981 GSM1385982 GSM1385983
do 
for j in `ls $i*`
do
awk '{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6"\t"$7}' $j > $j.bed &
done
done

# extract methylation profiles for cpgbody
for i in GSM1385973 GSM1385974 GSM1385975 GSM1385976 GSM1385977 GSM1385978 GSM1385979 GSM1385980 GSM1385981 GSM1385982 GSM1385983
do 
for j in `ls $i*.bed`
do
bedtools intersect -wa -a $j  -b cpgbody.bed > $j.cpg &
done
done

bedtools intersect -wa -a GSM1385983.dat.bed -b cpgbody.bed 

# decrease the sample size of the input to DSS 
for i in GSM1385973 GSM1385974 GSM1385975 GSM1385976 GSM1385977 GSM1385978 GSM1385979 GSM1385980 GSM1385981 GSM1385982 GSM1385983
do 
for j in `ls $i*.cpg`
do
awk '{print $1"\t"$2"\t$5"\t"$6}' $j > $j.dss &
done
done




for i in `seq 1..100000`
do
bedtools random -l 2 -n 596 -g mm9.len | bedtools intersect -wa -a - -b mm9.shelf.bed | wc -l
done



# get the raw data for the signficant site
source("http://bioconductor.org/biocLite.R")
library("heatmap.plus")
library("gplots")

setwd("/home/sguo/mice")
raw.sig<-read.table("matrix.detail.mice.alicey.txt",head=T,sep="\t")
data<-data.matrix(raw.sig[,c(7,4,8,9)])
clab<-c("x","129_mESC","B6_mESC","miPS_1","miPS_2","miPS_2","miPS_B3","SCNT_B12","SCNT_NB3","SCNT_P7C","SCNT_P8B")
colnames(data)=clab[c(7,4,8,9)]
pdf("mice.sig.pdf")
par(mar=c(1,1,1,1))
heatmap.plus(data,cex.main=0.2,cexRow=0.2,cexCol=0.9,col=redgreen(20));
dev.off()


'miPS_B3.BED.txt.trim',
'miPS_1E12P20.BED.txt.trim'
'SCNT_NB3.BED.txt.trim',
'SCNT_B12.BED.txt.trim'


setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\mice")



Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string="-wao"){
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

setwd("/home/shg047/monod/mhl")
d<-read.table("biomarker.24.gene.txt",head=F,sep="\t",as.is=T)[,1:4]

d1<-read.table("GEO.pvalue.dmr.rlt.total.txt",head=F,sep="\t",as.is=T)
d2<-read.table("coad.pvalue.dmr.rlt.total.txt",head=T,sep="\t",as.is=T)
d3<-read.table("luad.pvalue.dmr.rlt.total.txt",head=T,sep="\t",as.is=T)
d4<-read.table("lusc.pvalue.dmr.rlt.total.txt",head=T,sep="\t",as.is=T)
d5<-read.table("paad.pvalue.dmr.rlt.total.txt",head=T,sep="\t",as.is=T)


rlt1<-Rbedtools(functionstring="intersectBed",bed1=d,bed2=d1[,1:4],opt.string="-wao")

mat1<-d1[match(unique(rlt1[,8]),d1[,4]),]
mat2<-d2[match(unique(rlt1[,8]),d2[,1]),]
mat3<-d3[match(unique(rlt1[,8]),d3[,1]),]
mat4<-d4[match(unique(rlt1[,8]),d4[,1]),]
mat5<-d5[match(unique(rlt1[,8]),d5[,1]),]

rlt1<-na.omit(data.frame(mat1[,c(1:6,9,12)],mat2$pvalue,mat2$beta,mat3$pvalue,mat3$beta,mat4$pvalue,mat4$beta,mat5$pvalue,mat5$beta))
colnames(rlt1)<-c("CHR","Start","End","cgsite","Gene Symbol","Location","GEO.Pvalue","GEO.Delta","COAD.Pvalue","COAD.Delta","LUAD.Pvalue","LUAD.Detal","LUSC.Pvalue","LUSC.delta","PAAD.Pvalue","PAAD.delta")
write.table(rlt1,file="GEO.coad.luad.lusc.paad.beta.diff.txt",col.names=T,row.names=F,quote=F,sep="\t")



setwd("/home/sguo/monod/data")
data1<-read.table(file="GEO.coad.luad.lusc.paad.beta.diff.txt",head=T,as.is=T,sep="\t")

load("TCGA.COAD.450K.RData")
load("TCGA.LUAD.450K.RData")
load("TCGA.LUSC.450K.RData")
load("TCGA.PAAD.450K.RData")

v1<-match(data1[,4],rownames(coad))
v2<-match(data1[,4],rownames(luad))
v3<-match(data1[,4],rownames(lusc))
v4<-match(data1[,4],rownames(paad))

d1<-coad[v1,]
d2<-luad[v2,]
d3<-lusc[v3,]
d4<-paad[v4,]

function(data){
  u<-apply(data,1,function(x) sum(x[which(substr(colnames(d1),14,15)=="01")]>0.3))
  v<-apply(data,1,function(x) sum(x[which(substr(colnames(d1),14,15)=="11")]>0.3))
  
}




