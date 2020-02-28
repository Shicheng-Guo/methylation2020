setwd("/home/shg047/oasis/monod/HMHPlasma")
data=read.table("6-P-1.HMH.txt",head=T,sep="\t",check.names=F)
newdata=data[,4:ncol(data)]
newdata[newdata>0]=1

HeatMap<-function(data,cexRow = 0.2,cexCol = 0.9,Colv=T,Rowv=T){
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  colors <-greenred(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  ColSideColors=sidecol(colnames(data))
  heatmap.2(data,trace="none",cexRow = cexRow,cexCol = cexCol, ColSideColors=ColSideColors,density.info="none",col=colors,Colv=Colv,Rowv=Rowv,keysize=0.9, margins = c(5, 10))
}


newdata<-newdata[,-which(colnames(newdata)=="NC-P")]
newdata=newdata[order(colSums(newdata))]
# colon: 1351
# Lung: 1419 
# Almost Same. Chang their location
# which(colnames(newdata)=="Colon")
# which(colnames(newdata)=="Lung")
#Newdata=newdata
#Newdata[,9]=newdata[,8]
#Newdata[,8]=newdata[,9]
#colnames(Newdata)[9]=colnames(newdata)[8]
#colnames(Newdata)[8]=colnames(newdata)[9]

pdf("heatmap-no-NCP.pdf")
par(mfrow=c(3,5))
HeatMap(data.matrix(newdata),cexCol = 0.9,Colv=F)
dev.off()

pdf("heatmap-greenred-no-NCP.pdf")
HeatMap(data.matrix(newdata),cexCol = 0.9,Colv=F)
dev.off()

# Aim2: 
setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma")
pdf("heatmap-no-NCP-F15.pdf")
par(mfrow=c(3,5))
file=list.files(pattern="*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
data=read.table(file[i],head=T,sep="\t",check.names=F)
newdata=data[,4:ncol(data)]
newdata[newdata>0]=1
newdata<-newdata[,order(colSums(newdata))]
HeatMap(data.matrix(newdata),cexCol = 0.9,Colv=F)
print(i)
}
dev.off()


# Aim 2-1: heatmap for all the 15 samples.
setwd("/home/shg047/oasis/monod/HMHPlasma")
file=list.files(pattern="*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="pair-cancer-plasma.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)


# Aim2: heatmap for all the 15 samples.
setwd("/home/shg047/oasis/monod/HMHPlasma")
pdf("heatmap-no-NCP-F15.pdf")
par(mfrow=c(3,5))
file=list.files(pattern="*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[,-which(colnames(newdata)=="NC-P")]
  HeatMap(data.matrix(newdata),cexCol = 0.9,Colv=F)
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
dev.off()

write.table(yy,file="pair-cancer-plasma.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)
# newdata=newdata[order(colSums(newdata))]

# Aim3: Number of cancer specific methylation haplotype
setwd("/home/shg047/oasis/monod/HMHPlasma")
file=list.files(pattern="*.HMH.txt")
xse<-c()
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
 #  newdata<-newdata[,-c(which(colnames(newdata)=="NC-P"),which(colnames(newdata)=="WB"))]
  # newdata<-newdata[,-c(which(colnames(newdata)=="WB"))]
  # newdata<-newdata[,-c(which(colnames(newdata)=="WB"))]
  xsetmp<-which(rowSums(newdata)==2)
  newdata<-newdata[xsetmp,]
  xsetmpp<-paste(data[xsetmp,1],data[xsetmp,2],data[xsetmp,3],sep="-")
  xse<-c(xse,xsetmpp)
  print<-paste(file[i],length(which(rowSums(newdata)==2)),sep=" ")
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
# unique(xse)
xse251<-unique(xse)
length(xse251)
write.table(yy,file="tmp.txt",sep="\t",quote=F,col.names=NA,row.names=T)


# Aim3: Number of cancer specific methylation haplotype
setwd("/home/shg047/oasis/monod/HMHPlasma")
file=list.files(pattern="*.HMH.txt")
xse<-c()
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  xsetmp<-which(rowSums(newdata)==2)
  newdata<-newdata[xsetmp,]
  xsetmpp<-paste(data[xsetmp,1],data[xsetmp,2],data[xsetmp,3],sep="-")
  xse<-c(xse,xsetmpp)
  print<-paste(file[i],length(which(rowSums(newdata)==2)),sep=" ")
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
unique(xse)
xse251<-unique(xse)
length(xse251)
newdata<-newdata[,-which(colnames(newdata)=="NC-P")]

# Aim 4: Number of cancer specific methylation haplotype based non-pair data
setwd("/home/shg047/oasis/monod/HMHPlasma/Case3")
file=list.files(pattern="*.HMH.txt")
yy<-c()
xse<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  xsetmp<-which(rowSums(newdata)==1)
  xsetmpp<-paste(data[xsetmp,1],data[xsetmp,2],data[xsetmp,3],sep="-")
  xse<-c(xse,xsetmpp)
  xnum<-length(unique(data[xsetmp,1]))
  print<-paste(file[i],length(which(rowSums(newdata)==1)),xnum,sep=" ")
  print(print)
}
xse4244<-unique(xse)
length()
newdata<-newdata[,-which(colnames(newdata)=="NC-P")]

match(xse251,xse4244)

# Aim 5. Normal Plasma composition

setwd("/home/shg047/oasis/monod/HMHPlasma/NP-hapinfo")
file=list.files(pattern="*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file

setwd("/home/shg047/oasis/monod/HMHPlasma")
file=list.files(pattern="*.HMH.txt")
xse<-c()
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[,-c(which(colnames(newdata)=="NC-P"),which(colnames(newdata)=="WB"))]
  xsetmp<-which(rowSums(newdata)==2)
  newdata<-newdata[xsetmp,]
  xsetmpp<-paste(data[xsetmp,1],data[xsetmp,2],data[xsetmp,3],sep="-")
  xse<-c(xse,xsetmpp)
  print<-paste(file[i],length(which(rowSums(newdata)==2)),sep=" ")
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}

# Aim 5. Cancer Plasma composition

setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
file=list.files(pattern="*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)


x<-read.table("A.txt")
plot(x[,1],x[,2])
