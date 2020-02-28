setwd("/media/Home_Raid1/shg047/NAS3/HM450/TCGA/lihc")
target<-read.table("TCGA-LIHC.DMS.txt",sep="\t",head=T,check.names=F,as.is=T,row.names=18)
newtarget<-data.frame(CHR=target$V1,START=target$V2-250,END=target$V3+250,target)
pbmc<-read.table("/media/NAS3_volume2/shg047/HM450/TCGA/Normal.PBMC.GEO.HM450K.Beta.txt",sep="\t",row.names=1,head=T,check.names=F)
PBMC<-pbmc[match(rownames(newtarget),rownames(pbmc)),]
PBMCSubset<-subset(PBMC,mean<0.1 & median<0.1 & PBMC[,7]<0.15)
dim(PBMCSubset)

target<-data.frame(newtarget[match(rownames(PBMCSubset),rownames(newtarget)),],PBMCSubset,CPGSITE=rownames(PBMCSubset))
dim(target)
target<-subset(target,Statistic.pair>0.2)
dim(target)

write.table(target,file="TCGA-ESCA-PBMC.txt",col.names=T,row.names=F,sep="\t",quote=F)
bed<-data.frame(CHR=target[,1],START=target[,2]-250,END=target[,3]+250,ID=paste(target[,1],":",target[,2],"-",target[,3],sep=""))
write.table(bed,file="TCGA-ESCA-PBMC.bed",col.names=F,row.names=F,sep="\t",quote=F)



cd /media/Home_Raid1/shg047/NAS3/Roadmap/wig
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/NAS3/HM450/TCGA/lihc/TCGA-ESCA-PBMC.bed $i.out
echo $i
done


cd /media/Home_Raid1/shg047/NAS3/Roadmap/wig
for i in `ls *bw`
do
bigWigAverageOverBed $i /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed $i.out.mhb
echo $i
done


cd /media/Home_Raid1/shg047/NAS3/Roadmap/wig
for i in `ls *bw`
do
#bigWigAverageOverBed $i /media/Home_Raid1/shg047/NAS3/db/hg19/CpGI.hg19.bed4 $i.out.cpgi
#bigWigAverageOverBed $i /media/Home_Raid1/shg047/NAS3/db/hg19/CpG.Shore.hg19.bed4 $i.out.cpgshore
bigWigAverageOverBed $i /media/Home_Raid1/shg047/NAS3/db/hg19/Human-hg19-H3k27ac-PMID24119843.bed4 $i.out.h3k27ac
done


lial_Cells.Bisulfite-Seq.RM066.wig.gz.bw.out.h3k27ac

bigWigAverageOverBed GSM1127057_UCSF-UBC.Breast_Myoepithelial_Cells.Bisulfite-Seq.RM066.wig.gz.bw /media/Home_Raid1/shg047/NAS3/db/hg19/Human-hg19-H3k27ac-PMID24119843.bed4 GSM1127057_UCSF-UBC.Breast_Myoepithelial_Cells.Bisulfite-Seq.RM066.wig.gz.bw.out.h3k27ac



setwd("/media/Home_Raid1/shg047/NAS3/Roadmap/wig")
file=list.files(pattern="*.out.h3k27ac")

file=list.files(pattern="*.out.mhb")
data<-c()
for(i in file){
  tmp<-read.table(i,sep="\t",as.is=T,row.names=1)
  data<-cbind(tmp[,5],data)
}
rownames(data)<-rownames(tmp)
colName<-unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[2]))
Rlt<-apply(data,1,function(x) tapply(x,colName,mean))
colnames(Rlt)<-rownames(data)

RRlt<-t(apply(Rlt,2,quantile))
colnames(RRlt)<-paste("Normal",colnames(RRlt),sep="")

write.table(RRlt,file="/media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca/Target.Normal.Roadmap.meth.txt",row.names=T,col.names=NA,quote=F,sep="\t")
system("rm *.out")

fmax<-apply(Rlt,2,function(x) max(x)/sum(x))

Target<-data.frame(target,RRlt)
write.table(Target,file="/media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca/TCGA-ESCA-PBMC-NormalRoadmap.txt",row.names=F,col.names=T,quote=F,sep="\t")


# Enhancer Region()
for i in H3k4me1 H3k4me2 H3k4me3 H3k9ac H3k9me1 H3k9me3 H3k27ac H3k27me3 H3k36me3 H3k79me2 H4k20me1 
do
echo $i
bedtools intersect -wo -a TCGA-ESCA-PBMC.bed -b ~/NAS3/db/hg19/Human-hg19-$i-PMID24119843.bed | wc -l 
bedtools intersect -wo -a TCGA-ESCA-PBMC.bed -b ~/NAS3/db/hg19/Human-hg19-$i-PMID24119843.bed > TCGA-ESCA-PBMC.$i.bed
done

awk '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" TCGA-ESCA-PBMC-Histone.bed > TCGA-ESCA-PBMC-Histone.wig.bed
/media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca/TCGA-ESCA-PBMC-Histone.wig.bed

for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca/TCGA-ESCA-PBMC-Histone.wig.bed $i.out
done

? sort

target1<-read.table("/media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca/TCGA-ESCA-PBMC.bed",sep="\t",head=F,check.names=F,as.is=T)
target2<-read.table("/media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca/TCGA-ESCA.DMS.txt",sep="\t",head=T,check.names=F,as.is=T,row.names=18)

Target<-cbind(target1,RRlt)
colnames(Target)<-c("CHR","START","END",colnames(target2),colnames(RRlt))

write.table(Target,file="TCGA-ESCA-PBMC-Normal.txt",row.names=T,col.names=NA,quote=F,sep="\t")
