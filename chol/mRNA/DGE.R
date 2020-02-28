#######################################################################################################################
#######################################  CHOL DNA methylation Project ##############################################
#######################################################################################################################

DGEtest<-function(data,phen,compare.group=c("01","11")){
  output<-matrix(NA,dim(data)[1],5)   # set output matrix ()
  x1<-which(phen==compare.group[1])   # type 1, cancer or sensitive
  x2<-which(phen==compare.group[2])   # type 2, normal or resistant
  a.shapiro.pvalue<-c()
  b.shapiro.pvalue<-c()
  delta<-mean.A<-mean.B<-tmp.ttest<-tmp.wilcox<-c()
  Pvalue.ttest<-Pvalue.wilcox<-c()
  for(i in 1:nrow(data)){
    tmp.a.shapiro.pvalue<-try(shapiro.test(as.numeric(data[i,x1]))$p.value)
    tmp.b.shapiro.pvalue<-try(shapiro.test(as.numeric(data[i,x2]))$p.value)
    a.shapiro.pvalue<-c(a.shapiro.pvalue,tmp.a.shapiro.pvalue)
    b.shapiro.pvalue<-c(b.shapiro.pvalue,tmp.b.shapiro.pvalue)
    delta<-c(delta,mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2])))
    mean.A<-c(mean.A,mean(as.numeric(data[i,x1])))
    mean.B<-c(mean.B,mean(as.numeric(data[i,x2])))
    tmp.ttest<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2],paired=F)))
    Pvalue.ttest<-c(Pvalue.ttest,tmp.ttest$p.value)
    tmp.wilcox<-try(wilcox.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=F))
    Pvalue.wilcox<-c(Pvalue.wilcox,tmp.wilcox$p.value)
  }
  p.adjust.BH<-p.adjust(Pvalue.ttest,"BH")
  output<-data.frame(mean.A,mean.B,delta,Pvalue.ttest,p.adjust.BH,Pvalue.wilcox,a.shapiro.pvalue,b.shapiro.pvalue)
  rownames(output)<-rownames(data)
  return(output)
}


setwd("/media/Home_Raid1/shg047/work/tcga/chol/mrna")
file=list.files(pattern="*rsem.genes.normalized_results")
sampleinfo=read.table("unc.edu_CHOL.IlluminaHiSeq_RNASeqV2.1.0.0.sdrf.txt",head=T,sep="\t")
expDataSave1="/home/sguo/tcga/chol/TCGA.RNAseqV2.chol.RData"    # save RNAseqV2 matrix
exptxtSave1="/home/sguo/tcga/chol/TCGA.RNAseqV2.chol.Diff.txt"  # save signficant differential gene ouput

sampleid=as.character(sampleinfo[match(file,sampleinfo$Derived.Data.File),]$Comment..TCGA.Barcode.)
expdata<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=T,sep="\t",as.is=F)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  expdata<-cbind(expdata,tmp[,2])
  print(i)
}
rownames(expdata)=tmp[,1]
colnames(expdata)=sampleid
save(expdata,file="CHOL.mRNAseq.RData")
type<-substr(sampleid,14,15)
table(type)
expdata<-expdata+matrix(rnorm(length(expdata),0.0001,0.0001),dim(expdata)[1],dim(expdata)[2])   # row=gene, col=inv
rlt1<-DGEtest(expdata,phen=type,compare.group=c("01","11"))
write.table(rlt1,file="TCGA.CHOL.Differential.Gene.Expression.txt",col.names=T,row.names=T,quote=F,sep="\t")


# merge
setwd("/media/Home_Raid1/shg047/NAS1/tcga/chol/mrna")
rlt1<-read.table(file="../mrna/TCGA.CHOL.Differential.Gene.Expression.txt",row.names=1,sep="\t")
methgene<-read.table("methyGene.txt",head=T,sep="\t")
methgene<-read.table("../biomarker3.bed",head=F,sep="\t")
head(methgene)
output<-rlt1[match(methgene$V9,unlist(lapply(rownames(rlt1),function(x) unlist(strsplit(x,"[|]"))[1]))),]
xxrlt<-cbind(methgene,output)
write.table(xxrlt,file="Meth-Combin-TCGA.CHOL.Differential.Gene.Expression-time2.txt",col.names=NA,row.names=T,quote=F,sep="\t")
output2<-expdata[na.omit(match(methgene$V9,unlist(lapply(rownames(expdata),function(x) unlist(strsplit(x,"[|]"))[1])))),]

