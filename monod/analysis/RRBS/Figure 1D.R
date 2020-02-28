

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



# come to the 194 server
setwd("/home/sguo/monod/rrbs")
wgbsmb<-read.table("/home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",sep="\t")

##RRBS
load("RRBS.MethylationBlock.RData")
bed<-cor2bed(rownames(rrbsmb$hdrc))
rrbs<-data.frame(bed,rrbsmb$hdr,rrbsmb$hdrc)
rin<-Rbedtools(functionstring="intersectBed",rrbs,wgbsmb,opt.string="-wa -u")
rout<-Rbedtools(functionstring="intersectBed",rrbs,wgbsmb,opt.string="-wa -v")
rlt<-list()
rlt$rin<-rin
rlt$rout<-rout
save(rlt,file="RRBS.IN.OUT.Correlation.RData")  # download this file to desktop and plot the boxplot

##TCGA
setwd("/home/sguo/methylation")
load("COAD.mh.cor.RData")
wgbsmb<-read.table("/home/sguo/monod/phase3/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",sep="\t")

bed<-cor2bed(rownames(cor))
mh450<-data.frame(bed,cor)
rin<-Rbedtools(functionstring="intersectBed",mh450,wgbsmb,opt.string="-wa -u")
rout<-Rbedtools(functionstring="intersectBed",mh450,wgbsmb,opt.string="-wa -v")
rlt<-list()
rlt$rin<-rin
rlt$rout<-rout
save(rlt,file="mh450.IN.OUT.Correlation.RData")  # download this file to desktop and plot the boxplot

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\RRBS")
load("mh450.IN.OUT.Correlation.RData")
boxplot(rlt$rin[,4],rlt$rout[,4],outline=F,frame=T,names=c("IN","OUT"),col=c(2,4),lwd=2,cex.axis=1.4,ylab="Correlation",cex.lab=1.4)


head(rlt$rin)



