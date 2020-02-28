
library("Biostrings")
library("stringr")
library("randomForest")
library("impute")
library("rpart")
library("e1071")
library("biclust")

setwd("C:/Users/shicheng/Dropbox/Project/methylation/monod")

setwd("/home/shg047/monod/mhl")
sam<-read.table("../sampleinfo.sort.txt",sep="\t",as.is=T)
new<-cbind(sapply(sam[,1],function(x) unlist(strsplit(x,":"))[2]),sam[,2])
new<-data.frame(new)
write.table(newsam,file="../sampleinfo.sort.new.txt",sep="\t",quote=F,col.names=NA,row.names=T)

file1<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file2<-read.table("150209_BSPP_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file3<-read.table("140917_dRRBS_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file4<-read.table("141216_SeqCap_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file5<-read.table("150209_SeqCap_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)

cor1<-match(colnames(file1),new[,1])
cor2<-match(colnames(file2),new[,1])
cor3<-match(colnames(file3),new[,1])
cor4<-match(colnames(file4),new[,1])
cor5<-match(colnames(file5),new[,1])

lab1<-new[cor1,2]
lab2<-new[cor2,2]
lab3<-new[cor3,2]
lab4<-new[cor4,2]
lab5<-new[cor5,2]

table(new[cor1,2])
table(new[cor2,2])
table(new[cor3,2])
table(new[cor4,2])
table(new[cor5,2])

lab1<-lab1
file1<-file1
cor1<-cor1
dataset="RRBS"


rrbs<-pipeline(lab1=lab1,file1=file1,cor1=cor1,dataset="RRBS")
bspp<-pipeline(lab1=lab2,file1=file2,cor1=cor2,dataset="BSPP")
seqcap<-pipeline(lab1=lab5,file1=file5,cor1=cor5,dataset="SeqCap")

pipeline<-function(lab1=lab1,file1=file1,cor1=cor1,dataset="RRBS"){
rlt<-list()
u<-which(lab1=="Cancer")
v<-which(lab1=="Normal")
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
dim(outcomes)

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
best.mtry1=bestmtry$mtry1
rf.model4<-randomForest(t(newdata2),as.factor(lab1),mtry=best.mtry1,mtrees=500,importance=FALSE)
rlt$acc2<-senspeacc(rf.model4)

#importance2<-data.frame(ord=rownames(rf.model4$importance),importance=rf.model4$importance)
colnames(newdata)=colnames(file1)
Ofile5=paste(dataset,nrow(newdata1),"dataset","txt",sep=".")
Ofile6=paste(dataset,nrow(newdata2),"dataset","txt",sep=".")
write.table(newdata1,file=Ofile5,sep="\t",row.names=T,col.names=NA,quote=F)
write.table(newdata2,file=Ofile6,sep="\t",row.names=T,col.names=NA,quote=F)

return(rlt)
}


