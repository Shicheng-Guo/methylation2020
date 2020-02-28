# genome-wide DNA methylation comparsion between plasma and solid tissues.

dir1="/media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/";                   # RRBS1
dir2="/media/Ext12T/DD_Ext12T/MONOD/150209_SN216/BSPP/BAMfiles/";        # BSPP
dir3="/media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles/";        # RRBS2
dir4="/media/Ext12T/DD_Ext12T/MONOD/141216_HiSeqRapidRun/BAMfiles/";     # HiSeqRapidRun
dir5="/media/Ext12T/DD_Ext12T/MONOD/150209_SN216/SeqCap/BAMfiles/";      # SeqCap

setwd(dir1)
file1<-list.files(pattern="*.bam")
setwd(dir2)
file2<-list.files(pattern="*.bam")
setwd(dir3)
file3<-list.files(pattern="*.bam")
setwd(dir4)
file4<-list.files(pattern="*.bam")
setwd(dir5)
file5<-list.files(pattern="*.bam")

sam1<-unlist(unique(lapply(file1,function(x) strsplit(x,"[._]")[[1]][1])))
sam2<-unlist(unique(lapply(file2,function(x) strsplit(x,"[._]")[[1]][1])))
sam3<-unlist(unique(lapply(file3,function(x) strsplit(x,"[._]")[[1]][1])))
sam4<-unlist(unique(lapply(file4,function(x) strsplit(x,"[._]")[[1]][1])))
sam5<-unlist(unique(lapply(file5,function(x) strsplit(x,"[._]")[[1]][1])))

#the barcode system is not consistant, for example 6-P-1 and 6P-1. Remove all the hyphen. 
sam11<-gsub("-","",sam1)
sam21<-gsub("-","",sam2)
sam31<-gsub("-","",sam3)
sam41<-gsub("-","",sam4)
sam51<-gsub("-","",sam5)

length(sam11)
length(sam21)
length(sam31)
length(sam41)
length(sam51)

sam<-sam21[na.omit(match(sam11,sam21))]
sam<-sam21[na.omit(match(sam51,sam21))]
write.table(sam,file="shared.sample.rrbs.bspp.txt",col.names=NA,row.names=NA,sep="\t",quote=F)


substr(sam11,1,2)
sam<-unique(c(sam11,sam21,sam31,sam41,sam51))
NC<-sam[grep("NC",sam)]
NC<-unique(gsub("CP","C",NC))
length(NC)

samfilt<-function(x){
  rlt1<-x[grep("6T",x)]
  rlt2<-x[grep("6P",x)]
  rlt3<-x[grep("7T",x)]
  rlt4<-x[grep("7P",x)]
  rlt5<-x[grep("PCT",x)]
  rlt6<-x[grep("PCP",x)]
  rlt<-data.frame(rlt1,rlt2,rlt3,rlt4,rlt5,rlt6)
  rlt
}

filt1<-samfilt(sam11)
filt2<-samfilt(sam21)
filt3<-samfilt(sam31)
filt4<-samfilt(sam41)
filt5<-samfilt(sam51)

# only dir1 and dir3, aka RRBS and dRRBS have samples both for solid tissues and corresponding plasma. 
setwd("/home/shg047/monod/mhl")
rrbs<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,row.names=1,sep="\t",as.is=T,check.names=F)
rrbs[23,11]
rrbs[23,12]
# here, we can find the mhl is not stable in different run for the same sample

colnames(rrbs)
rownames(rrbs)
# extract paired samples
pairsam<-function(x){
  rlt<-list()
  sam<-unlist(unique(lapply(x,function(x) strsplit(x,"[._]")[[1]][1])))
  t<-sam[grep("-T-",sam)]
  t2p<-gsub("T","P",sam[grep("-T-",sam)])
  p<-sam[grep("-P-",sam)]
  pair<-c(t,p[match(t2p,p)])
  sam2<-unlist((lapply(x,function(x) strsplit(x,"[._]")[[1]][1])))
  size<-length(match(pair,sam2))/2
  rlt$position<-match(pair,sam2)[unlist(lapply(1:size,function(x) c(x,x+size)))]
  rlt$sampleid<-sam2[rlt$position]
  rlt
}

rrbs_pair<-rrbs[,pairsam(colnames(rrbs))$position]
colnames(rrbs_pair)<-pairsam(colnames(rrbs))$sampleid

par(mfrow=c(3,5))

cor<-c()
for(i in 1:15){
  tmp<-cor(rrbs_pair[,2*i-1],rrbs_pair[,2*i])
  cor<-c(cor,tmp)
}


################################################################################################
###################Correlation analysis between tissues and plasma##############################
################################################################################################
setwd("C:/Users/shicheng/Dropbox/Project/methylation/monod")
setwd("/home/shg047/monod/mhl")
load("rrbs_pair.RData")
jpeg("cor.3.title.distribution.jpeg")
plot(cor,type="l",col=rainbow(10)[1],lwd=2,xlab="Sample ID",ylab="Correlation",ylim=c(0.1,0.7))
ci95(cor)
z<-1
for(j in seq(0.2,0.9,by=0.2)){
z<-z+1
row<-c()
for(i in 1:30){
  tmp<-which(rrbs_pair[,i]>j)
  row<-c(row,tmp)
  row<-unique(row)
}
cor1<-c()
for(i in 1:15){
  tmp<-cor(rrbs_pair[row,2*i-1],rrbs_pair[row,2*i])
  cor1<-c(cor1,tmp)
}
lines(cor1,col=rainbow(10)[2*z],lwd=2)
ci95(cor1)
}
legend("bottomright",legend=c("0.0","0.2","0.4","0.6","0.8"),col=rainbow(10)[c(1,4,6,8,10)],lty=1,lwd=2,bty="n")
dev.off()


ci95<-function(x){
error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
m<-round(mean(x),2)
d<-round(mean(x)-error,2)
u<-round(mean(x)+error,2)
rlt<-paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
print(rlt)
}


sum(rrbs_pair[,1]>0.1)/nrow(rrbs_pair)
sum(rrbs_pair[,2]>0.1)/nrow(rrbs_pair)
sum(rrbs_pair[,1]==0)/nrow(rrbs_pair)
sum(rrbs_pair[,2]==0)/nrow(rrbs_pair)


# co-methylation regions both in tissues and plasma simuteniously. 
load("rrbs_pair.RData")
jpeg("overlap.number.jpeg")
plot(1,1,type="n",col=rainbow(10)[1],lwd=2,xlab="Sample ID",ylab="Ratio",ylim=c(0.1,0.58),xlim=c(1,15))
ratio<-c()
z<-0
  for(j in seq(0,0.9,0.1)){
    z<-z+1
    u<-c()
  for(i in 1:15){
    x1<-which(rrbs_pair[,2*i-1]>j)
    x2<-which(rrbs_pair[,2*i]>j)
    r<-sum(x1 %in% x2)/length(unique(c(x1)))
    u<-c(u,r)
  }
  lines(u,col=rainbow(10)[2*z],lwd=2)
  ratio<-rbind(ratio,u)
  }
legend("bottomright",legend=c("0.0","0.2","0.4","0.6","0.8"),col=rainbow(10)[c(2,4,6,8,10)],lty=1,lwd=2,bty="n")
dev.off()


load("rrbs_pair.RData")
for(i in 1:15){
  loci<-which(rrbs_pair[,2*i-1]*rrbs_pair[,2*i]>0)
  decrease<-which(rrbs_pair[loci,2*i-1]-rrbs_pair[loci,2*i]>0.2)
  increase<-which(rrbs_pair[loci,2*i-1]-rrbs_pair[loci,2*i]< -0.2)
  length(decrease)
  length(increase)
  par(mfrow=c(1,2))
  hist(rrbs_pair[,2*i]-rrbs_pair[,2*i-1],breaks=300,xlab="Delta(Tissue - Plasma)",main="",cex=2)
  hist(rrbs_pair[loci,2*i]-rrbs_pair[loci,2*i-1],breaks=300,xlab="Delta(Tissue - Plasma)",main="")
  lines(density(rrbs_pair[loci,2*i]-rrbs_pair[loci,2*i-1]))
}

rrbs and colon cancer
5358 regions or fragment showed methylated status both in tissues and plasma. 


################################################################################################
######################################### Random and selective test#############################
################################################################################################
setwd("/home/shg047/monod/mhl")
load("rrbs_pair.RData")
region<-c()
for(i in 1:15){
  for(j in 1:nrow(rrbs_pair)){
    if(rrbs_pair[j,2*i-1]*rrbs_pair[j,2*i]>0){
    region<-c(region,j)
    }
  }
}
region<-unique(region)
newdata<-rrbs_pair[region,]
sum<-0
for(i in 1:15){
  for(j in 1:nrow(newdata)){
    if(newdata[j,2*i-1]*newdata[j,2*i]>0){
      sum<-sum+1
    }
  }
}
p=sum/(nrow(newdata)*ncol(newdata))
  Pvalue<-c()
  regions<-c()
  for(i in 1:nrow(newdata)){
    region<-0
    for(j in 1:15){
      if(newdata[i,2*j-1]*newdata[i,2*j]>0){
        region<-region+1
      }
    }
    ptmp<-binom.test(region,15,p,alternative="two.sided")$p.value
    Pvalue<-c(Pvalue,ptmp)
    regions<-c(regions,region)
}

sig<-which(Pvalue<0.05/length(Pvalue))
newdata2<-cbind(newdata[sig,],Pvalue[sig])

head(newdata2)

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

sigbed<-cor2bed(rownames(newdata2))
sigbed<-data.frame(sigbed,Pvalue[sig])
sigene<-Rbedtools(bed1=sigbed,bed2=refgene,opt="-wao")
write.table(sigene,file="over-slected.tissue.plasma.gene.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(newsigene,file="new.over-slected.tissue.plasma.gene.txt",col.names=F,row.names=F,quote=F,sep="\t")
head(refgene)




observed <- sort(Pvalue)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
pdf("qqplot.pdf", width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

