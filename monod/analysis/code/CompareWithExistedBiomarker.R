

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


setwd("/home/shg047/oasis/monod/mhl")
data<-read.table("colon.cancer.plasma.vs.normal.plasma.significant.txt",head=T)
bed<-cor2bed(rownames(data))
newdata<-cbind(bed,data)
write.table(newdata,"lung.cancer.plasma.vs.normal.plasma.significant.txt",col.names=F,row.names=F,sep="\t",quote=F)
system("bedtools intersect -wao -a lung.cancer.plasma.vs.normal.plasma.significant.txt -b /home/shg047/oasis/db/hg19.refGene.bed > lung.cancer.sig.mhl.sorted.Annotated.txt")


load("MONOD-Apr6.MHL.RData")
sum(is.na(data))/(nrow(data)*ncol(data))
# re-group/collect the samples #colon cancer
Group<-colnames(data)
# Group[grep("CTR|NC-|Pregn",Group)]<-"NP"
Group[grep("NC-",Group)]<-"NP-Kun"
Group[grep("CTR|Pregn",Group)]<-"NP-Dennis"
Group[grep("6-P-",Group)]<-"CCP"
Group[grep("6-T-|HCT116|Colon_Tumor_Primary|CTT-|metastasis_colon|tumor_colon",Group)]<-"CCT"
Group[grep("normal_colon|N37-Colon|SG-01",Group)]<-"NCT"

Group[grep("7-P-",Group)]<-"LCP"
Group[grep("7-T-|adenocarcinoma_lung|tumor_lung",Group)]<-"LCT"
Group[grep("LG-01|N37-Lung|normal_lung",Group)]<-"NLT"

Group[grep("PC-P-",Group)]<-"PCP"
Group[grep("PC-T-",Group)]<-"PCT"
Group[grep("PA-01|N37-Pancreas",Group)]<-"NPT"

Group[grep("WB-",Group)]<-"WB"
Group[grep("STL|N37-|methylC-",Group)]<-"ONT"
newdata<-data
colnames(newdata)<-Group
newdata<-data.frame(newdata,check.names=F)

row<-grep("chr1:25258448-25258484",rownames(newdata))
col<-c(grep("NP-Dennis",colnames(newdata)),grep("NP-Kun",colnames(newdata)),grep("CCP",colnames(newdata)))
col
newdata[row,col]


biomark<-read.table("colonCancerBiomarker.txt")
sig<-read.table("lung.cancer.sig.mhl.sorted.Annotated.txt",head=T,sep="\t")
match(sig[,13],biomark[,1])









