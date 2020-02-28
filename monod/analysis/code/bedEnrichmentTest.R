## usage: 
setwd("/home/sguo/monod/tcga")
rlt<-bedEnrichmentTest(bed1="Conservative.methylation.block.0.6.bed.txt",bed2="WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",iteration=1000000)
save(rlt,file="tcga.wgbs.overlap.RData")
##


bedEnrichmentTest<-function(bed1,bed2,assembly="hg19",iteration=1000000){
  rlt<-list()
  hsa.chr<-paste("chr",c(1:22,"X","Y"),sep="")
  hg19.size<-c("249250621","243199373","198022430","191154276","180915260","171115067",
               "159138663","146364022","141213431","135534747","135006516","133851895",
               "115169878","107349540","102531392","90354753","81195210","78077248",
               "59128983","63025520","48129895","51304566","155270560","59373566")
  hg38.size<-c("248956422","242193529","198295559","190214555","181538259","170805979",
               "159345973","145138636","138394717","133797422","135086622","133275309",
               "114364328","107043718","101991189","90338345","83257441","80373285",
               "58617616","64444167","46709983","50818468","156040895","57227415")
  mm.chr<-paste("chr",c(1:19,"X","Y"),sep="")
  mm9.size<-c("197195432","181748087","159599783","155630120","152537259","149517037",
              "152524553","131738871","124076172","129993255","121843856","121257530",
              "120284312","125194864","103494974","98319150","95272651","90772031",
              "61342430","166650296","15902555")
  mm10.size<-c("195471971","182113224","160039680","156508116","151834684","149736546",
               "145441459","129401213","124595110","130694993","122082543","120129022","120421639",
               "124902244","104043685","98207768","94987271","90702639","61431566","171031299","91744698")
  hg19.chrom.sizes<-data.frame(hsa.chr,hg19.size)
  hg38.chrom.sizes<-data.frame(hsa.chr,hg38.size)
  mm9.chrom.sizes<-data.frame(mm.chr,mm9.size)
  mm10.chrom.sizes<-data.frame(mm.chr,mm10.size)
  write.table(hg19.chrom.sizes,file="hg19.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(hg38.chrom.sizes,file="hg38.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(mm9.chrom.sizes,file="mm9.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(mm10.chrom.sizes,file="mm10.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  
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
  
  command=paste("bedtools intersect -wa -u","-a",bed1,"-b",bed2,">","output1",sep=" ")
  try(system(command))
  count1<-read.table("output1")
  observe<-nrow(count1)
  try(system(paste("rm","output1",sep=" " )))
  try(system(paste("rm","output2",sep=" " )))
  chrom.size=paste(assembly,"chrom.sizes",sep=".")
  for (i in 1:iteration){
    command=paste("bedtools shuffle","-i",bed1,"-g",chrom.size,"| bedtools intersect -wa -u -a -","-b",bed2,"|wc -l >> output2",sep=" ")
    try(system(command))
    print(i)
  }
  count2<-read.table("output2")
  
  print(paste("The expected overlap region number is",mean(count2[,1]),sep=" "))
  print(paste("The observed overlap region number is",observe,sep=" "))
  
  pvalue<-1-sum(count2[,1] < observe)/iteration
  if(pvalue==0){
    pvalue=paste("<",1/iteration,sep=" ");
  }else{
    pvalue=paste("=",pvalue,sep=" ")
  }
  print(paste("The P-value of the target bed region overlapped with source bed, Pvalue",pvalue,sep=""))
  try(system("mv output2 expected.overlab.counts.number"))
  
  pdf("enrichment.hist.box.pdf")
  hist(count2[,1],breaks=25,col="red",main="",
       ylab="Counts")
  boxplot(count2[,1], at=15,horizontal=TRUE,
          outline=F,boxwex=50, frame=F,border=3,col = "green1", add = TRUE,lwd=6)
  dev.off()
  
  rlt$expectation<-count2[,1]
  rlt$observation<-observe
  rlt$pvalue=pvalue
  return(rlt)
}

