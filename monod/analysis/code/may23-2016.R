

############################
# Composition of HMH in cancer plamsa for Colon cancer
setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
# file <- Sys.glob("6-P*.HMH.txt")
file=paste("6-P-",1:5,".HMH.txt",sep="")
zz<-c()
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
zz<-rbind(zz,yy)
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Colon")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Intestine")]==0),] 
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)
yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="6-T")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Colon")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Intestine")]==0),]  
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)
write.table(zz,file="cancdre.plasma.pairwise-colon-2.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)

#############################################
# Composition of HMH in cancer plamsa for Lung cancer
setwd("/oasis/tscc/scratch/shg047/monod/HMHPlasma/Case4")
# file <- Sys.glob("6-P*.HMH.txt")
file=paste("7-P-",1:5,".HMH.txt",sep="")
zz<-c()

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
zz<-rbind(zz,yy)

yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)


yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Lung")]==0),]
  
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)


yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="7-T")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Lung")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)

write.table(zz,file="cancdre.plasma.pairwise-lung-2.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)




#############################################
# Composition of HMH in cancer plamsa for Normal Plasma
setwd("/home/shg047/oasis/monod/HMHPlasma/NP-hapinfo")
file <- sort(Sys.glob("NC-P-*HMH.*"))
# file=paste("7-P-",1:5,".HMH.txt",sep="")
zz<-c()

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
zz<-cbind(zz,yy)

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
zz<-cbind(zz,yy)


yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Lung")]==0),]
  
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)


yy<-c()
for(i in 1:length(file)){
  data=read.table(file[i],head=T,sep="\t",check.names=F)
  newdata=data[,4:ncol(data)]
  newdata[newdata>0]=1
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="WB")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="7-T")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="NC-P")]==0),]
  newdata<-newdata[which(newdata[,which(colnames(newdata)=="Lung")]==0),]
  y<-colSums(newdata)
  yy<-rbind(yy,y)
  print(i)
}
rownames(yy)<-file
zz<-rbind(zz,yy)

write.table(zz,file="normal.plasma.componment.txt",sep="\t",quote=F,col.names=NA,row.names=T)






