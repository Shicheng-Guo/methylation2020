

args <- commandArgs(trailingOnly = TRUE)
file1<-args[1]
file2<-args[2]

setwd("/home/sguo/monod/data")
library("Biostrings")
library("stringr")
library("randomForest")
library("impute")
library("rpart")
library("e1071")
library("biclust")

RINfun=function(yorig){
  yranks=rank(yorig)
  tempp=(yranks-.5)/(length(yranks))
  return(qnorm(tempp))
}
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
cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}
pipeline<-function(lab1=lab1,file1=file1,cor1=cor1,dataset="RRBS"){
  rlt<-list()
  u<-which(lab1=="Cancer")
  v<-which(lab1=="Normal")
  rlt$cancer=u
  rlt$normal=v
  rlt$rawdata=file1
  keep<-c()
  for(i in 1:nrow(file1)){
    p1<-sum(file1[i,u]>0)
    p2<-sum(file1[i,v]==0)
    print(i)
    if(p2==length(v) & p1>0){
      tmp<-c(i,p1)
      keep<-rbind(keep,tmp)
    }  
  }
  rlt$totalnum<-nrow(file1)
  rlt$hyperfreq<-keep
  rlt$hypernum<-nrow(keep)
  
  # hypermethylation frequency in cancer samples for the loci which is totally un-methylation in normal samples. 
  Ofile1<-paste(dataset,"hypermethylation.freq.hist","jpeg",sep=".")
  jpeg(Ofile1)
  hist(keep[,2]/length(u),main="",xlab="Hypermethylation Frequency in Cancer",lwd=1.5,cex.lab=1.5,cex.axis=1.5)
  dev.off()
  
  # function used to get best mtry and mtree
  data<-as.matrix(file1[keep[,1],])
  bestmtry<-tuneRF(x=t(data), y=as.factor(lab1), ntreeTry=1000, stepFactor=2, improve=0.05,trace=TRUE, plot=TRUE, doBest=T)
  best.mtry=bestmtry$mtry
  senspeacc<-function(x){
    acc=(x$confusion[1,1]+x$confusion[2,2])/sum(x$confusion)
    sen=(x$confusion[1,1])/(x$confusion[1,1]+x$confusion[1,2])
    spe=(x$confusion[2,2])/(x$confusion[2,1]+x$confusion[2,2])
    c(sen,spe,acc)
  }
  err.rate<-c()
  senspacc<-c()
  for(mtrees in seq(250,5000,by=250)){
    rf.model1<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=mtrees)
    tmp1<-mean(rf.model1$err.rate[,1])
    tmp2<-senspeacc(rf.model1)
    err.rate<-c(err.rate,tmp1)
    senspacc<-c(senspacc,tmp2)
  }
  bestmtrees<-seq(250,5000,by=250)[which.min(err.rate)]
  
  # do the random forest based on above best mtry and mtree
  # informative loci which is useful in Random Forest model based on methylated regions in cancer tissues. 
  rf.model<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=bestmtrees,importance=FALSE)
  importance<-data.frame(ord=rownames(rf.model$importance),importance=rf.model$importance)
  num<-sum(importance[,2]>0)
  rlt$acc1<-senspeacc(rf.model)
  
  
  newdata1<-data[order(importance[,2],decreasing=T)[1:num],]
  keep1<-c()
  for(i in 1:nrow(newdata1)){
    p1<-sum(newdata1[i,u]>0)
    p2<-sum(newdata1[i,v]==0)
    print(i)
    if(p2==length(v) & p1>0){
      tmp<-c(i,p1)
      keep1<-rbind(keep1,tmp)
    }  
  }
  rlt$inforfreq<-keep1
  rlt$infornum<-nrow(keep1)
  
  Ofile2<-paste(dataset,"RF.informative.loci.freq.hist.in.cancer","jpeg",sep=".")
  jpeg(Ofile2)
  hist(keep1[,2]/length(u),main="",xlab="Hypermethylation Frequency in Cancer",lwd=1.5,cex.lab=1.5,cex.axis=1.5,col="red")
  dev.off()
  
  Ofile3<-paste(dataset,"informative.loci.importance.distribution.jpeg",sep=".")
  jpeg(Ofile3)
  plot(importance[order(importance[,2],decreasing=T)[1:num],2],main="",xlab="importance")
  dev.off()
  
  for(j in 1:num){
    a<-sum((importance[order(importance[,2],decreasing=T)[1:j],2])^2)
    b<-sum((importance[order(importance[,2],decreasing=T)[1:num],2])^2)
    if(a/b>0.8){
      bestnum<-j
      break
    } 
  }
  
  newdata2<-data[order(importance[,2],decreasing=T)[1:bestnum],]
  keep2<-c()
  for(i in 1:nrow(newdata2)){
    p1<-sum(newdata2[i,u]>0)
    p2<-sum(newdata2[i,v]==0)
    print(i)
    if(p2==length(v) & p1>0){
      tmp<-c(i,p1)
      keep2<-rbind(keep2,tmp)
    }  
  }
  rlt$top80freq<-keep2
  Ofile4<-paste(dataset,"RF.top80.informative.loci.freq.hist.in.cancer","jpeg",sep=".")
  jpeg(Ofile4)
  hist(keep2[,2]/length(u),main="",xlab="Hypermethylation Frequency in Cancer",lwd=1.5,cex.lab=1.5,cex.axis=1.5)
  dev.off()
  
  bestmtry1<-tuneRF(x=t(newdata2), y=as.factor(lab1), ntreeTry=5, stepFactor=2, improve=0.05,trace=TRUE, plot=TRUE, doBest=T)
  best.mtry1=bestmtry1$mtry
  rf.model4<-randomForest(t(newdata2),as.factor(lab1),mtry=best.mtry1,mtrees=500,importance=FALSE)
  rlt$acc2<-senspeacc(rf.model4)
  
  #importance2<-data.frame(ord=rownames(rf.model4$importance),importance=rf.model4$importance)
  colnames(newdata1)=colnames(file1)
  colnames(newdata2)=colnames(file1)
  Ofile5=paste(dataset,nrow(newdata1),"dataset","txt",sep=".")
  Ofile6=paste(dataset,nrow(newdata2),"dataset","txt",sep=".")
  write.table(newdata1,file=Ofile5,sep="\t",row.names=T,col.names=NA,quote=F)
  write.table(newdata2,file=Ofile6,sep="\t",row.names=T,col.names=NA,quote=F)
  
  return(rlt)
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
Rbedtoolsort<-function(functionstring="bedtools sort",bed1,opt.string=""){
  #create temp files
  a.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-i",a.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}
DMRtest<-function(case,control,outputfilename="output.result"){
  rlt<-list()
  data<-data.matrix(data.frame(case,control))
  x1<-ncol(case)
  x2<-ncol(control)
  x3<-ncol(data)
  # wilcox test and median, beta, fold fold 
  pvalue<-apply(data,1,function(x) p<-wilcox.test(x[1:x1],x[(x1+1):x3],na.omit=T)$p.value)
  med1<-apply(data,1,function(x) median(x[1:x1],na.rm=T))
  med2<-apply(data,1,function(x) median(x[(x1+1):x3],na.rm=T))
  beta<-med1-med2
  ratio<-med1/med2
  hypercase<-apply(data,1,function(x) sum(x[1:x1]>0.3,na.rm=T)/length(na.omit(x[1:x1])))
  hypercon<-apply(data,1,function(x) sum(x[(x1+1):x3]>0.3,na.rm=T)/length(na.omit(x[(x1+1):x3])))
  output<-data.frame(pvalue,med1,med2,beta,ratio,hypercase,hypercon)
  sig<-subset(output,abs(beta)>0.1 & ratio>1.2 & hypercase>0.4)
  output1=paste(outputfilename,".total.txt",sep="")
  output2=paste(outputfilename,".sig.txt",sep="")
  write.table(output,file=output1,col.names=NA,row.names=T,quote=F,sep="\t")
  write.table(sig,file=output2,col.names=NA,row.names=T,quote=F,sep="\t")
  return(output)
  rlt$pvalue=output
  rlt$sig=sig
}

load(file1)
tcgaDMRtest(data=coad,fileout="coad")