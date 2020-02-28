bed2cor<-function(bed){
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],"-",unlist(strsplit(x,"\t"))[3],sep="")})
  return(cor)
}

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

setwd("/home/sguo/monod/data/geo")
f<-read.table("/home/sguo/annotation/GPL13534.sort.bed")
f<-subset(f,f[6]=="Island")
f1<-read.table("share.hdrc.GSE35069.GSE41169.GSE42861.txt",sep="\t")
f2<-read.table("conservative.methylation.block.C11.0.6.bed.txt",sep="\t")

m1<-Rbedtools(functionstring="intersectBed",f1,f,opt.string="-wa -u")
m2<-Rbedtools(functionstring="intersectBed",f2,f,opt.string="-wa -u")
mm1<-Rbedtools(functionstring="intersectBed",m1,m2,opt.string="-wa -u")

dim(f1)
dim(f2)
dim(m1)
dim(m2)

m3<-Rbedtools(functionstring="intersectBed",f1,f2,opt.string="-wa -u")
m4<-Rbedtools(functionstring="intersectBed",m3,f,opt.string="-wa -u")
dim(m3)
dim(m4)

# random sampling PBMC and TCGA methylation block and estimate the expectation of the overlap
file<-list.files(pattern="*Share.High.HDRC.0.6.bed.txt")
t1<-read.table(file[1],as.is=T)
t2<-read.table(file[2],as.is=T)
t3<-read.table(file[3],as.is=T)

tmp1<-unique(rbind(t1,t2,t3)) # 984 
tmp2<-unique(read.table("/home/sguo/methylation/conservative.methylation.block.0.6.bed"))
dim(tmp1)
dim(tmp2)

b<-c()
for(i in 1:10000){
b1<-tmp1[sample(1:nrow(tmp1),52),]
b2<-tmp2[sample(1:nrow(tmp2),72),]  
b3<-try(Rbedtools(functionstring="intersectBed",b1,b2,opt.string="-wa -u"))
b<-c(b,nrow(b3))
}
write.table(b,file="background.20150818.txt",col.names=F,row.names=F,quote=F)  # save to wiki


library("VennDiagram")
# Venn for PBMC
area1=138
area2=590
area3=909
n12=123
n13=121
n23=528
n123=119
pdf("a.pdf")
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, cex=3,cat.cex=2,col.lab="white",cat.pos=c(350,0,180),category =c("GSE35069","GSE41169","GSE42861"),col=c(2:4),fill=2:4,lwd=2,ind = T, list.order = 1:3)
dev.off()
? draw.triple.venn
# Venn for PBMC and TCGA
area2=155
area1=119
cross.area=51
draw.pairwise.venn(area1, area2, cross.area, category = c("PBMC","TCGA-Normal"),col=2:3,cat.pos=c(330,0),fill=2:3,cat.cex=2.5,cex=3.5)

area1=155
area2=276
cross.area=55
draw.pairwise.venn(area2, area1, cross.area, category = c("PBMC","TCGA-Cancer"),rotation.degree = 180,col=2:3,cat.pos=c(180,180),fill=2:3,cat.cex=2.5,cex=3.5)

# Venn for PBMC, TCGA-Normal, TCGA-Cancer
area1=155
area2=276
area3=119
n12=113  #nrow(rlt3)
n13=55   #nrow(rlt2)
n23=51   #nrow(rlt1)
n123=43  #nrow(rlt4)
pdf("a.pdf")
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, cex=3,cat.cex=1.5,col.lab="white",cat.pos=c(350,0,180),category =c("TCGA-Normal","TCGA-Cancer","PBMC"),col=c(2:4),fill=2:4,lwd=2,ind = T, list.order = 1:3)
dev.off()


# CpG island
x1<-c(67,33,84)
33/67
33/84
# non-CpG island
x2<-c(52,18,71)
18/52
18/71

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod")
load("simulation.20150818-1.RData")
par(mar=c(3,4,1,1))
hist(b,col="red",main="",xlab="",ylim=c(0,4000))






