
sam<-read.table("../sampleinfo.sort.txt",sep="\t",as.is=T)
new<-cbind(sapply(sam[,1],function(x) unlist(strsplit(x,":"))[2]),sam[,2])
new<-data.frame(new)
refgene<-read.table("../../annotation/hg19.bed")


bedwithgap<-function(bed,gap){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
}


file1<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file2<-read.table("150209_BSPP_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file3<-read.table("140917_dRRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file4<-read.table("141216_SeqCap_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file5<-read.table("150209_SeqCap_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)

cor1<-match(colnames(file1),new[,1])
cor2<-match(colnames(file2),new[,1])
cor3<-match(colnames(file3),new[,1])
cor4<-match(colnames(file4),new[,1])
cor5<-match(colnames(file5),new[,1])

lab1<-new[cor1,2]
lab2<-new[cor2,2]
lab3<-new[cor3,2]
lab4<-new[cor4,2]
lab5<-new[cor5,2]

load(file="rrbs.rlt.RData")
load(file="bspp.rlt.RData")
load(file="seqcap.rlt.RData")




hybed1<-data.frame(cor2bed(rownames(rrbs$rawdata)[rrbs$hyperfreq[,1]]),rrbs$hyperfreq[,2])
hybed2<-data.frame(cor2bed(rownames(seqcap$rawdata)[seqcap$hyperfreq[,1]]),seqcap$hyperfreq[,2])
hybed3<-data.frame(cor2bed(rownames(bspp$rawdata)[bspp$hyperfreq[,1]]),bspp$hyperfreq[,2])


avg<-c()
for(gap in seq(1,300,by=1)){
  gap=25
  hybed1wgap<-bedwithgap(hybed1,gap)
  hybed2wgap<-bedwithgap(hybed2,gap)
  hybed3wgap<-bedwithgap(hybed3,gap)
  head(rownames(rrbs$rawdata)[rrbs$hyperfreq[,1]])
  head(hybed1)
  head(hybed1wgap)
  int1<-Rbedtools(bed1=hybed1wgap,bed2=hybed2wgap,opt.string="-wo")
  nrow(int1)
  int2<-Rbedtools(bed1=int1,bed2=hybed3wgap,opt.string="-wo")
  nrow(int2)
  bio<-paste(int2[,10],":",int2[,11]+gap,"-",int2[,12]-gap,sep="")
  bsppDiag<-t(bspp$rawdata[match(bio,rownames(bspp$rawdata)),])
  u<-which(lab2=="Cancer")
  tmp<-c(gap,sum(rowSums(bsppDiag[u,])>0)/length(u))
  avg<-rbind(avg,tmp)
  print(tmp)
}