


setwd("/home/sguo/mice/geo/trim")
phen<-read.table("phen.txt")
load("MethylationMatrixHuamnReprogramming.RData")
file=list.files(pattern="*dat.RData")
filename<-unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[1]))
colnames(data)<-phen[match(filename,phen[,1]),2]
rowname<-paste("chr",rownames(data),sep="")
rownames(data)<-rowname



input<-read.table("")

for i in `ls *fixed.bed.txt`
do
awk '{print $1":"$2"\t"$4}' $i > $i.sim &
  done


awk '{print $1":"$2"\t"$4}' GSM1385973.dat > GSM1385973.dat.fixed.bed.sim



file=list.files(pattern="*fixed.bed.sim.RData")
load(file[1])
shareloc<-tmp[,1]
for(i in 2:length(file)){
  load(file[i])
  shareloc<-shareloc[which(shareloc %in% tmp[,1])]
  print(file[i])
}
save(shareloc,file="shareloc.RData")

setwd("/home/sguo/mice/geo/trim")
load("shareloc.RData")
file=list.files(pattern="*fixed.bed.sim.RData")
matrix<-c()
for(i in 1:length(file)){
  load(file[i])
  matrix<-cbind(matrix,tmp[match(shareloc,tmp[,1]),2])
  print(file[i])
}
rlt<-data.frame(shareloc,matrix)
colnames(rlt)<-c("LOC",file)
save(rlt,file="MethylationMatrixHuamnReprogramming.RData")









load("GSM1385973.dat.RData")
load("GSM1385981.dat.fixed.bed.sim.RData")

s1<-paste(dat1[,1],dat1[,2],sep=":")
s2<-dat1[,5]/dat1[,6]
tmp<-data.frame(s1,s2)
save(tmp,file="GSM1385973.dat.fixed.bed.sim")



setwd("/home/sguo/annotation/chipseq")



d1<-read.table("output_DMS.tsv",head=T)
d2<-read.table("DMS.human.txt")

dd1<-paste(d1[,1],d1[,2],sep=":")


dd2<-paste(d2[,1],d2[,2],sep=":")

d2[which(dd2 %in% dd1 ==T),]

awk '{print "chr"$1"\t"$2"\t"$2+1}' DMS.human.txt > DMS.human.txt.bed

bedtools interesect -wa -a /home/sguo/annotation/chipseq/wgEncodeRegTfbsClusteredWithCellsV3.bed -b /home/sguo/mice/geo/trim/DMS.human.txt.bed 



chr<-paste("chr",sapply(strsplit(as.character(rlt[,1]),":"),function(x) unlist(x)[1]),sep="")
start<-as.numeric(sapply(strsplit(as.character(rlt[,1]),":"),function(x) unlist(x)[2]))
end<-start+1




data<-data.frame(chr,start,end,rlt[,2:ncol(rlt)])
write.table(data,file="human.reprogramming.methylation.data.bed",sep="\t",row.names=T,col.names=NA,quote=F)

phen<-read.table("phen.txt")
colnames(raw)<-colnames(data)

sam<-sapply(strsplit(colnames(raw)[4:ncol(raw)],"[.]"),function(x) unlist(x)[1])

phen[match(sam,phen[,1]),2]



chr<-paste("chr",sapply(strsplit(as.character(shareloc),":"),function(x) unlist(x)[1]),sep="")
start<-as.numeric(sapply(strsplit(as.character(shareloc),":"),function(x) unlist(x)[2]))
end<-start+1

data<-data.frame(chr,start,end,matrix)
colnames(data)<-c("chr","start","end",file)
sam<-sapply(strsplit(colnames(data)[4:ncol(data)],"[.]"),function(x) unlist(x)[1])
colnames(data)<-c("chr","start","end",sam)
write.table(data,file="human.reprogramming.methylation.data.bed",sep="\t",row.names=F,col.names=F,quote=F)


bedtools intersect -wa -a human.reprogramming.methylation.data.newbed -b Sig1.bed > sig1.raw.single.data.txt

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\mice")
datak<-read.table("sig1.raw.single.data.txt",sep="")
colnames(datak)<-colnames(data)
x<-datak[,2]-26169470
y1<-datak[,5]
y2<-rowMeans(datak[,c(6,7)])
y3<-rowMeans(datak[,c(4,8:10)])
y4<-rowMeans(datak[,c(11:14)])
newdata<-data.frame(x,y1,y2,y3,y4)

pdf("methylation.profile.2.pdf")
plot(0,0,xlim=c(26169470,max(datak[,2])),ylim=c(0,1),type="n",xlab="Genome Postion",ylab="methylation level")
smoothingSpline = smooth.spline(datak[,2], newdata[,2], spar=1.4)
lines(smoothingSpline,col=2-1,lty=2,lwd=4)
for(i in 3:ncol(newdata)){
  smoothingSpline = smooth.spline(datak[,2], newdata[,i], spar=1.4)
  lines(smoothingSpline,col=i-1,type="l",pch=2,lwd=4)
  }
legend("bottomleft",c("Fibroblast","HESO","SCNT","iPS"),col=1:4,lty=1,bg="transparent",bty="n",lwd=4,cex=0.75)
dev.off()





