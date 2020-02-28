

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



setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\layer_specfic_mhl")
/home/shg047/monod/dec/Table.GSI.layer.mhl.WGBS.Remove.H1.WBC.rlt.txt
data=read.table("/home/shg047/monod/dec/Table.GSI.layer.mhl.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
head(data)
newdata=subset(data,GSI>0.6)
table(subset[,2])
mesoderm<-subset(newdata,group=="Mesoderm")
ectoderm<-subset(newdata,group=="Ectoderm")
endoderm<-subset(newdata,group=="Endoderm")

mesodermBed<-cor2bed(mesoderm[,1])
ectodermBed<-cor2bed(ectoderm[,1])
endodermBed<-cor2bed(endoderm[,1])


# change to laptop to plot histgram
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\layer_specfic_mhl")
data=read.table("Table.GSI.layer.mhl.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
newdata=subset(data,GSI>0.6)
write.table(newdata,file="Supplementary Table. layer specific MHL regions.txt",sep="\t",row.names=F)

pdf("hist.layer.specfic.pdf")
hist(data[,3],breaks=30,col="green",ylim=c(0,4000),xlab="Layer Specfic Index",main="")
dev.off()


Rbedtools(functionstring="intersectBed",bed1=mesodermBed,bed2=ectodermBed,opt.string="-wa -u")


table<-c(99,66,9)
names(table)<-c("Ectoderm","Endoderm","Mesoderm")
pdf("layer.specfic.number.pdf",width=5,height=5)
barplot(table,ylim=c(0,100),col="green")
dev.off()
