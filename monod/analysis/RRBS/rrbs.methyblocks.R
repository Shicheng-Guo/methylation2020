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
cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}
bed2cor<-function(bed){
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],"-",unlist(strsplit(x,"\t"))[3],sep="")})
  cor<-gsub("[ ]","",cor)
  return(cor)
}





################################################ Methylation Haplotype Block for RRBS dataset ##############################

setwd("/home/sguo/monod/rrbs")
load("encode.rrbs.pos.freq.RData")
dim("encode.rrbs.pos.freq")  # 2646999
load("Encode.RRBS.map.RData")
dim(map) # 866979
load("pos.RData")
upos<-unique(pos)
bed<-cor2bed(upos)
map<-data.frame(bed,upos)
write.table(map,file="rrbs.encode.bedmap",sep="\t",quote=F,col.names=T,row.names=F) 
# sort the bedmap file

map<-read.table("rrbs.encode.bedmap",head=T,sep="\t",as.is=T)
data<-read.table("rrbs.encode.data.txt",sep="\t",as.is=T,row.names=1,head=T)
rrbsmb<-mbsearch(data,map)
save(rrbsmb,file="RRBS.MethylationBlock.RData")

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\RRBS")


load("RRBS.MethylationBlock.RData")
rlt<-hist(rrbsmb$hdr[,3],breaks=50,xlim=c(0,60),main="")
plot(rlt$breaks[1:10],rlt$counts[1:10],type="o",lwd=3,cex=1.25,xlab="Number of CpG",ylab="Frequency")

pdf("a.pdf",width=2,height=2)
hist(rrbsmb$hdrc[,1],breaks=100,main="",xlab="Pearson correlation")
dev.off()


setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\RRBS")

setwd("G:\\monod\\encode")
load("RRBS.MethylationBlock.RData")
bed<-cor2bed(rownames(rrbsmb$hdrc))
dim(rrbsmb$hdr)
dim(rrbsmb$hdrc)

summary(lm(rrbsmb$hdrc[,1]~as.numeric(as.character(rrbsmb$hdr[,3]))))
summary(lm(rrbsmb$hdrc[,1]~as.numeric(as.character(rrbsmb$hdr[,4]))))
summary(lm(rrbsmb$hdrc[,1]~as.numeric(as.character(rrbsmb$hdr[,5]))))
summary(lm(rrbsmb$hdrc[,1]~as.numeric(as.character(rrbsmb$hdr[,6]))))

sum(rrbsmb$hdrc[,1]>0.6,na.rm=T)/length(rrbsmb$hdrc[,1]) # 

rrbs<-data.frame(rrbsmb$hdr,rrbsmb$hdrc)
head(rrbs)
library("ggplot2")
qplot(CpGNumber,ratio,data=rrbs)
qplot(RegionLength,ratio,data=rrbs)
qplot(CpGRatio,ratio,data=rrbs)

pdf("adc.pdf")
qplot(ratio,RegionLength,data=rrbs)
qplot(ratio,CpGNumber,data=rrbs)
dev.off()
getwd()


subset6<-cor2bed(rownames(subset(rrbs,ratio>0.6)))
subset7<-cor2bed(rownames(subset(rrbs,ratio>0.7)))
subset8<-cor2bed(rownames(subset(rrbs,ratio>0.8)))
subset9<-cor2bed(rownames(subset(rrbs,ratio>0.9)))
subset95<-cor2bed(rownames(subset(rrbs,ratio>0.95)))
subset99<-cor2bed(rownames(subset(rrbs,ratio>0.99)))


dim(subset6)  # 23517
dim(subset7)  # 18071
dim(subset8)  # 10986
dim(subset9)  # 3395
dim(subset95)  # 672
dim(subset99)  # 25

write.table(subset6,file="RRBS.MHB.0.6.matrix.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(subset7,file="RRBS.MHB.0.7.matrix.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(subset8,file="RRBS.MHB.0.8.matrix.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(subset9,file="RRBS.MHB.0.9.matrix.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(subset95,file="RRBS.MHB.0.95.matrix.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(subset99,file="RRBS.MHB.0.99.matrix.txt",sep="\t",quote=F,col.names=F,row.names=F)


bedtools intersect -wa -a RRBS.MHB.0.6.matrix.txt -b /home/shg047/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l   # 1856
bedtools intersect -wa -a RRBS.MHB.0.7.matrix.txt -b /home/shg047/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l   # 1703
bedtools intersect -wa -a RRBS.MHB.0.8.matrix.txt -b /home/shg047/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l   # 1259
bedtools intersect -wa -a RRBS.MHB.0.9.matrix.txt -b /home/shg047/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l   # 476
bedtools intersect -wa -a RRBS.MHB.0.95.matrix.txt -b /home/shg047/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l   # 106
bedtools intersect -wa -a RRBS.MHB.0.99.matrix.txt -b /home/shg047/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l   # 1





head(subset1)
save(subset1,file="Encode.RRBS.MHB.RData") # 23517

qplot(ratio,RegionLength,data=subset1)

shapiro.test(subset1$RegionLength)
ksnormTest(subset1$RegionLength)
getwd()

rlt<-data.frame(bed,rrbsmb$hdrc,rrbsmb$hdr)
write.table(rlt,file="rrbs.")

rrbs<-na.omit(rrbs)
smoothScatter(x=rrbs$ratio,y=rrbs$RegionLength,nbin=128,ylim=c(0,600))

? smoothScatter

library("fBasics")

wgbsmb<-read.table("/home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",sep="\t")

bedtools intersect -wa -a -b /home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed



cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

bed<-cor2bed(as.character(map[,4]))

map<-read.table("Encode.RRBS.map.sort",as.is=T,sep="\t")
dist<-c()
for(i in paste("chr",c(1:22,"X","Y"),sep="")){
tmp<-subset(map,map[,1]==i)
for(j in 2:nrow(tmp)){
n<-tmp[j,3]-tmp[j-1,3] 
dist<-c(dist,n)
}
print(paste(i,"is on the calculation..."),sep=" ")
}


# level of methylation block identified in TCGA dataset and Encode dataset overlapped with GWBS methylation blocks.
setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\monod\\RRBS")
hist1<-read.table("Conservative.methylation.block.0.6.bed.txt.number")
hist(hist1[,1],breaks=25,col="red",main="",ylab="Counts")
boxplot(hist1[,1], at=15,horizontal=TRUE,outline=F,boxwex=50, frame=F,border=3,col = "green1", add = TRUE,lwd=6)



setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\RRBS")
rrbs<-read.table("rrbs.output2.num")
hm450<-read.table("tcga.output2.num")

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\mhb")
pdf("mh450.random.pdf")
boxplot(hm450[,1],ylim=c(1,60),cex=0.8,yaxt="n",lwd=2,outline=T,frame=F,width=1,xlab="Array")
axis(2,at=seq(0,60,by=10),labels=c(0,50,100,150,1000,1200,1400),cex.axis=1.25)
text(x=1,y=50,labels="*",cex=3)
dev.off()

pdf("RRBS-encode.random.pdf")
boxplot(rrbs[,1],ylim=c(30,320),cex=0.8,yaxt="n",lwd=2,outline=T,frame=F,width=1,xlab="RRBS")
axis(2,at=seq(30,320,by=50),labels=c(150,300,450,8000,9000,10000),cex.axis=1.25)
text(x=1,y=231,labels="*",cex=3)
dev.off()


# methylation block for RRBS dataset
mbsearch<-function(data,sortedmapfile,window=100,mincpgnumber=4){
  # rowname of data is cpg site/identifier and column name is sample id
  # map file is bed file and the fourth column is cpg/identifier
  library("impute")
  output<-list()
  data<-RawNARemove(data)
  data<-impute.knn(data.matrix(data))$data
  map<-sortedmapfile
  map<-map[map[,4] %in% rownames(data),]
  newdata<-data[match(map[,4],rownames(data)),]
  if(nrow(map)==nrow(newdata)){
    a<-map[,2]
    i=1
    tmp<-c()
    rlt<-c()
    rowname<-c()
    index<-0
    while(i < length(a)){
      start=a[i]
      end=a[i+1]
      if(end-start<window && end-start>0){
        tmp<-c(tmp,i)
        i=i+1
      }else{
        if(length(tmp)>mincpgnumber){
          index=index+1
          tmp<-c(min(tmp),max(tmp),length(tmp),a[max(tmp)]-a[min(tmp)],round(length(tmp)/(a[max(tmp)]-a[min(tmp)]),4),round((a[max(tmp)]-a[min(tmp)])/(length(tmp))))
          rlt<-rbind(rlt,tmp)
          tmp2<-paste(map[tmp[1],1],":",map[tmp[1],2],"-",map[tmp[2],2],sep="")
          rowname<-c(rowname,tmp2)
        }
        tmp<-c()
        i=i+1
      }
    }
    rownames(rlt)<-rowname
    colnames(rlt)<-c("rowstart","rowstart","CpGNumber","RegionLength","CpGRatio","AverageGap")
    cor<-c()
    for(j in 1:nrow(rlt)){
      cor1<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],]),use="complete.obs")) # cancer
      cor<-c(cor,cor1)
    }
    cor<-data.frame(cor)
    rownames(cor)<-rowname
    colnames(cor)<-c("ratio")
    
    output$hdr<-rlt
    output$hdrc<-cor
    
    return(output)
  }
}
