

library("impute")

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
} 
bedwithgap<-function(bed,gap){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
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

bed2cor<-function(bed){
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],"-",unlist(strsplit(x,"\t"))[3],sep="")})
  cor<-gsub("[ ]","",cor)
  return(cor)
}



f1<-read.table(file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
# heat only have 3 high GSI regions
sum(table(subset(f1,GSI>=0.6)[,2]))
f1<-subset(f1,GSI>0.6)
bed1<-cor2bed(f1$region)

f2<-read.table("RRBS_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(f2)
f2<-RawNARemove(f2,missratio=0.4)
f2<-impute.knn(data.matrix(f2))$data
s2<-unlist(lapply(colnames(f2),Name2Group))
colnames(f2)<-s2
remove=c(s2[grep("-T-",s2)])
f2<-f2[,-match(remove,s2)]
bed2<-cor2bed(rownames(f2))
bed21<-Rbedtools(functionstring="intersectBed",bed2,bed1,opt.string="-wa -u")
cor2<-bed2cor(bed21)
f2<-f2[match(cor2,rownames(f2)),]
# Stage I: random forest with normal plasma

x<-t(f2)
y<-as.factor(sapply(colnames(f2),function(x) substr(x,1,2)))
library("randomForest")
fit<-randomForest(scale(x),y,importance=T)
top1<-match(names(sort(fit$importance[,6],decreasing=T)[1:ncol(x)]),colnames(x))
fit1<-randomForest(scale(x)[,top1],y,importance=T)
fit1

for(j in seq(5,150,by=5)){
top2<-match(names(sort(fit$importance[,6],decreasing=T)[1:50]),colnames(x))
fit2<-randomForest(scale(x)[,top2],y,importance=T)
print(c(j,fit2$confusion[,5]))
}

# Stage II: random forest with normal plasma
f2=f2[,-grep("NC",colnames(f2))]
x<-t(f2)
y<-as.factor(sapply(colnames(f2),function(x) substr(x,1,2)))
library("randomForest")
fit<-randomForest(scale(x),y,importance=T)
for(j in seq(5,100,by=5)){
  top2<-match(names(sort(fit$importance[,5],decreasing=T)[1:j]),colnames(x))
  fit2<-randomForest(scale(x)[,top2],y,importance=T)
  print(c(j,fit2$confusion[,4]))
}

top2<-match(names(sort(fit$importance[,5],decreasing=T)[1:60]),colnames(x))
fit2<-randomForest(scale(x)[,top2],y,importance=T)

### Seq-Cap
f3<-read.table("/home/shg047/monod/oct/data/WGBS_SeqCap_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
f3<-RawNARemove(f3,missratio=0.35)
f3<-impute.knn(data.matrix(f3))$data
bed3<-cor2bed(rownames(f3))
bed31<-Rbedtools(functionstring="intersectBed",bed3,bed1,opt.string="-wa -u")
cor3<-bed2cor(bed31)
f3<-f3[match(cor3,rownames(f3)),]

s3<-unlist(lapply(colnames(f3),Name2Group))
colnames(f3)<-s3

x<-t(f3)
y<-as.factor(sapply(colnames(f3),function(x) substr(x,1,2)))
library("randomForest")
fit<-randomForest(scale(x),y,importance=T)
for(j in seq(2,30,by=1)){
  top2<-match(names(sort(fit$importance[,6],decreasing=T)[1:j]),colnames(x))
  fit2<-randomForest(scale(x)[,top2],y,importance=T)
  print(c(j,fit2$confusion[,5]))
}

top2<-match(names(sort(fit$importance[,5],decreasing=T)[1:13]),colnames(x))
fit2<-randomForest(scale(x)[,top2],y,importance=T)


f3=f3[,-grep("NC",colnames(f3))]
x<-t(f3)
y<-as.factor(sapply(colnames(f3),function(x) substr(x,1,2)))
library("randomForest")
fit<-randomForest(scale(x),y,importance=T)
for(j in seq(2,30,by=1)){
  top2<-match(names(sort(fit$importance[,5],decreasing=T)[1:j]),colnames(x))
  fit2<-randomForest(scale(x)[,top2],y,importance=T)
  print(c(j,fit2$confusion[,4]))
}

top2<-match(names(sort(fit$importance[,5],decreasing=T)[1:3]),colnames(x))
fit2<-randomForest(scale(x)[,top2],y,importance=T)




Name2Group<-function(Name1){
  Name1=unlist(strsplit(Name1,"[.]"))[1]
  Name1=unlist(strsplit(Name1,"_"))[1]     
  Name1<-gsub("6-P","CC-P",Name1)
  Name1<-gsub("6P","CC-P",Name1)
  Name1<-gsub("7P","LC-P",Name1)
  Name1<-gsub("PCP","PC-P",Name1)
  Name1<-gsub("7-P","LC-P",Name1)
  Name1<-gsub("6-T","CC-T",Name1)
  Name1<-gsub("7-T","LC-T",Name1)
  Name1<-gsub("-frozen","-f",Name1)
  Name1<-gsub("-FFPE","-F",Name1)
  Name1<-gsub("-100ng","-2U",Name1)
  Name1<-gsub("-5ng","-1U",Name1)
  Name1<-gsub("CTT","CC-T",Name1)
  return(Name1)
}












