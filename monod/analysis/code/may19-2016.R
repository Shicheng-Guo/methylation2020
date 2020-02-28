
# Aim 6-1. Composition of HMH in cancer plamsa. 

setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
# file <- Sys.glob("6-P*.HMH.txt")
file=paste("6-P-",1:5,".HMH.txt",sep="")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
#  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
#  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.pairwise-1.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)


setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
# file <- Sys.glob("6-P*.HMH.txt")
file=paste("6-P-",1:5,".HMH.txt",sep="")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.pairwise-2.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)


setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
# file <- Sys.glob("6-P*.HMH.txt")
file=paste("6-P-",1:5,".HMH.txt",sep="")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Colon")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Intestine")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.pairwise-3.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)


setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
# file <- Sys.glob("6-P*.HMH.txt")
file=paste("6-P-",1:5,".HMH.txt",sep="")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="6-T")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Colon")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Intestine")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.pairwise-4.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)



# Aim 5-2. HMH not in Normal Plasma

setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
file <- Sys.glob("*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.nonp-2.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)

setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
file <- Sys.glob("*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.nonp-3.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)


setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
file <- Sys.glob("*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
write.table(yy,file="cancdre.plasma.nonp-3.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)




# Aim 5-2. HMH not in any normal tissues
# mv 6-T, 7-T, 8PC-T to bak
setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
file <- Sys.glob("6-*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1 
  newdata<-newdata[,-c(3:4)]
  newdata<-subset(newdata,newdata[,1]==1 & newdata[,2]==1)
  newdata<-newdata[which(rowSums(newdata)==2),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
head(yy)
write.table(yy,file="cancdre.plasma.nonp-nont-6.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)


file <- Sys.glob("7-*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1  
  newdata<-newdata[,-c(1,4)]
  newdata<-subset(newdata,newdata[,1]==1 & newdata[,2]==1)
  newdata<-newdata[which(rowSums(newdata)==2),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
head(yy)
write.table(yy,file="cancdre.plasma.nonp-nont-7.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)

file <- Sys.glob("PC-*.HMH.txt")
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1  
  newdata<-newdata[,-c(1,2)]
  newdata<-subset(newdata,newdata[,1]==1 & newdata[,11]==1)
  newdata<-newdata[which(rowSums(newdata)==2),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
yy
write.table(yy,file="cancdre.plasma.nonp-nont-PC.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)







