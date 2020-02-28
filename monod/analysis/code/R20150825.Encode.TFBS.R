

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

bedEnrichmentTest<-function(bed1,bed2,assembly="hg19",iteration=1000){
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
  excl.file="CRGmapability.hg19.exclude.bed";
  for (i in 1:iteration){
    command=paste("bedtools shuffle","-i",bed1,"-g",chrom.size,"-excl",excl.file,"| bedtools intersect -wa -u -a -","-b",bed2,"|wc -l >> output2",sep=" ")
    try(system(command))
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
  
  rlt$expectation<-count2[,1]
  rlt$observation<-observe
  rlt$pvalue=pvalue
  return(rlt)
}


setwd("/home/sguo/annotation/hg19/")
file<-list.files(pattern="encode.*bed")
bed<-read.table("/home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",sep="\t",as.is=T)

count<-c()
for(i in 1:length(file)){
x<-read.table(file[i],sep="\t",as.is=T)
y<-Rbedtools(functionstring="intersectBed",bed1=x,bed2=bed,opt.string="-wa -u")
r<-nrow(y)/nrow(x)
tf<-unlist(strsplit(file[i],"[.]"))[2]
tmp<-c(tf,nrow(x),nrow(y),r)
count<-rbind(count,tmp)
}

pvalue.enrich<-c()
for(i in 1:length(file)){
rlt<-bedEnrichmentTest(bed1=file[i],bed2="/home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",iteration=100)
tmp<-c(file[i],rlt$pvalue)
pvalue.enrich<-rbind(pvalue.enrich,tmp)
print(tmp)
}

for i in seq 1,10,100
do
bedtools shuffle -i encode.ATF2.hg19.bed -g hg19.chrom.size -excl CRGmapability.hg19.exclude.bed | bedtools intersect -wa -u -a - -b /home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l
done
bedtools intersect -wa -u -a encode.ATF2.hg19.bed -b /home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l


pdf("barplot.pdf")
barplot(as.numeric(count[order(count[,4]),4]),col =heat.colors(nrow(count)),ylim=c(0,0.1))
dev.off()
write.table(count[order(count[,4]),],file="counts.txt",sep="\t",quote=F)


x<-c(0.888,0.48,0.27,0.28,0.26,0.26,0.26,0.076,0.046,0.0418,0.0474,0.036,0.054,0.031,0.071,0.04,0.052,0.036)
names(x)<-c("LDA","LOCKS","H1 DMR","ME DMR","NPC DMR","MSC DMR","TBL DMR","CpGI","H1 UMR","H1 LMR","ME UMR","ME LMR","NPC UMR","NPC LMR","MSC UMR","MSC LMR","TBL UMR","TBL LMR")
barplot(x, col =rainbow(12),horiz=F,cex.names = 0.7,ylim=c(0,1))
