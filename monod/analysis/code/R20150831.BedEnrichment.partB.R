###################################################################################
###   Title: Bed enrichment analysis to comprehensive reference database. (part B)
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 9/31/2015
###################################################################################


setwd("/home/sguo/monod/tcga")
source("bedEnrichmentTest.R")

setwd("/home/sguo/annotation/hg19/")
ref=list.files(pattern="*hg19.bed")
# ref=c("LAD.hg19.bed","LOCK.hg19.bed","Coloncancer.small.DMR.PMC3145050.hg19.bed","VMR.charm.PMID20844285.hg19.bed",
#   "Hic.common.boundary.hESC.IMR90.hg19.bed","Hic.topological.domain.IMR90.hg19.bed","Hic.topological.domain.hESC.hg19.bed",
#   "CpGI.hg19.bed","CpG.Shore.hg19.bed")

#### Part B (bed1 <-bed2, bed2<-bed1)
setwd("/home/sguo/monod/phase3/partB")
observation<-c()
expectation<-c()
ratio<-c()
fc<-c()
pvalue<-c()

for(i in ref){
  bed2="WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed"
  bed1=paste("/home/sguo/annotation/hg19/",i,sep="")
  excl="/home/sguo/annotation/hg19/CRGmapability.hg19.exclude.bed"
  assembly="hg19"
  nbed1=read.table(bed1)
  nbed2=read.table(bed2)
  rlt<-bedEnrichmentTest(bed1=bed1,bed2=bed2,excl=excl,iteration=1000,assembly)
  observation=c(observation,rlt$observation)
  expectation=c(expectation,mean(rlt$expectation))
  fc=c(fc,rlt$observation/mean(rlt$expectation))
  ratio=c(ratio,rlt$observation/nrow(nbed1))
  pvalue=c(pvalue,rlt$pvalue)
}
out<-data.frame(ref,observation,expectation,ratio,fc,pvalue)
write.table(out,file="WGBS_pooled_mappable_bins.bedenrichment.result.partB.txt",sep="\t",quote=F,row.names=F,col.names=T)


