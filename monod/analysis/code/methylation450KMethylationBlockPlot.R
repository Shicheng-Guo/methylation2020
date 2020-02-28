setwd("/home/sguo/methylation")
RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
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

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

library("impute")
map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")
file=list.files(pattern="*pair.RData")
load(file[j])
cancer<-substr(file[j],1,4)
data<-RawNARemove(data)
data<-impute.knn(data)$data
map<-map[map[,4] %in% rownames(data),]
newdata<-data[match(map[,4],rownames(data)),]

chr1:2345333-2345410

load("COAD.mh.cor.RData")
load("COAD.mh.rlt.RData")

j<-match("chr1:17445929-17446087",rownames(rlt))

cor1<-cor(t(newdata[rlt[j,1]:rlt[j,2],seq(1,ncol(newdata),by=2)]),use="complete.obs") # cancer
cor2<-cor(t(newdata[rlt[j,1]:rlt[j,2],seq(2,ncol(newdata),by=2)]),use="complete.obs") # normal

x1<-as.vector(newdata[rlt[j,1]:rlt[j,2],seq(1,ncol(newdata),by=2)])
x2<-as.vector(newdata[rlt[j,1]:rlt[j,2],seq(2,ncol(newdata),by=2)])

colnames(x)<-c("Cancer","Normal")
pdf("heatmap.b.pdf")
boxplot(x[,1],x[,2],outline=F,xlab="",cex.axis=1.25,names=c("Cancer","Normal"),col=c(2,3))
dev.off()

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod")
load("X.RData")

? boxplot

pdf("heatmap.b.pdf")
heatmap(newdata[rlt[j,1]:rlt[j,2],seq(1,ncol(newdata),by=2)],key=T)
heatmap(newdata[rlt[j,1]:rlt[j,2],seq(2,ncol(newdata),by=2)],key=T)
dev.off()



M1<-rltt$cor1
M2<-rltt$cor2
M1[lower.tri(M1)] <- NA
M2[lower.tri(M2)] <- NA
col=colorRampPalette(c("white", "red"))(20) 
pdf("figure.pdf")
image(M1,col =col,frame=F,xaxt="n",yaxt="n")
image(M2,col =col,frame=F,xaxt="n",yaxt="n")
dev.off()

library("grDevices")
plot(1:20,pch=20,col=col)
col=colorRampPalette(c("white", "red"))(20) 
M <- matrix(runif(100),10,10)
M[lower.tri(M)] <- NA
image(M,col = col,frame=F,xaxt="n",yaxt="n")
