library("IRanges")
library("Biostrings")
library("stringr")
library("randomForest")
library("impute")
library("rpart")
library("e1071")
library("biclust")


setwd("C:/Users/shicheng/Dropbox/Project/methylation/monod")
sam<-read.table("sampleinfo.sort.txt",sep="\t",as.is=T)
new<-cbind(sapply(sam[,1],function(x) unlist(strsplit(x,":"))[2]),sam[,2])
new<-data.frame(new)
write.table(new,file="../sampleinfo.sort.new.txt",sep="\t",quote=F,col.names=NA,row.names=T)

file<-list.files(pattern="*stringent_mhl_matrix_*")[1:6]
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

u<-which(lab1=="Cancer")
v<-which(lab1=="Normal")

keepnrow<-function(){
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
keep
}

hist(keep[,2]/length(u),xlab="Hypermethylation Frequency in Cancers",ymain="",cex=1.5)

data<-as.matrix(file1[keep[,1],])
dim(data)

bestmtry<-tuneRF(x=t(data), y=as.factor(lab1), ntreeTry=500, stepFactor=2, improve=0.05,trace=TRUE, plot=TRUE, doBest=T)
best.mtry=bestmtry$mtry

? tuneRF

rf.model1<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=250)
rf.model2<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=500)
rf.model3<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=750)
rf.model4<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=1000)
rf.model5<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=3000)
rf.model6<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=5000)
rf.model7<-randomForest(t(data),as.factor(lab1),mtry=best.mtry,mtrees=7000)
mean(rf.model1$err.rate[,1])
mean(rf.model2$err.rate[,1])
mean(rf.model3$err.rate[,1])
mean(rf.model4$err.rate[,1])
mean(rf.model5$err.rate[,1])
mean(rf.model6$err.rate[,1])
mean(rf.model7$err.rate[,1])

rf.model1
rf.model2
rf.model3
rf.model4
rf.model5
rf.model5
rf.model7

senspacc<-function(x){
  acc=(x$confusion[1,1]+x$confusion[2,2])/sum(x$confusion)
  sen=(x$confusion[1,1])/(x$confusion[1,1]+x$confusion[1,2])
  spe=(x$confusion[2,2])/(x$confusion[2,1]+x$confusion[2,2])
  c(sen,spe,acc)
}

senspacc(rf.model1)
senspacc(rf.model2)
senspacc(rf.model3)
senspacc(rf.model4)
senspacc(rf.model5)
senspacc(rf.model6)
senspacc(rf.model7)

rf.model4<-randomForest(t(data),as.factor(lab1),mtry=bestmtry$mtry,mtrees=1000,importance=FALSE)
importance<-data.frame(ord=rownames(rf.model4$importance),importance=rf.model4$importance)
num<-sum(importance[,2]>0)
num
jpeg("Capseq.importance.distribution.jpeg")
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
newdata<-data[order(importance[,2],decreasing=T)[1:bestnum],]

bestmtry2<-tuneRF(x=t(newdata), y=as.factor(lab1), ntreeTry=5, stepFactor=2, improve=0.05,trace=TRUE, plot=TRUE, doBest=T)
rf.model4<-randomForest(t(newdata),as.factor(lab1),mtry=bestmtry2$mtry,mtrees=1000,importance=FALSE)
importance2<-data.frame(ord=rownames(rf.model4$importance),importance=rf.model4$importance)
colnames(newdata)=colnames(file1)

write.table(newdata,file="Capseq.516.newdata.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(newdata,file="Capseq.21.newdata.txt",sep="\t",row.names=T,col.names=NA,quote=F)

write.table(newdata,file="RRBS.1592.newdata.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(newdata,file="RRBS.289.newdata.txt",sep="\t",row.names=T,col.names=NA,quote=F)


methylationFreq<-function(x){
keep<-c()
for(i in 1:nrow(newdata)){
  p1<-sum(newdata[i,u]>0)
  p2<-sum(newdata[i,v]==0)
  print(i)
  if(p2==8 & p1>0){
    tmp<-c(i,p1)
    keep<-rbind(keep,tmp)
  }  
}
keep
}


mean<-mean(keep[,2]/31)
sd<-sd(keep[,2]/31)
IQR<-IQR(keep[,2]/31)


wilcox.test(x[u],x[v])$p.value
p.adjust<-p.adjust(p, "BH") 

pdf("keep.proportion.pdf")
hist(keep[,2],breaks=20,main="",xlab="Cumulation of hypermethylated samples",xlim=c(0,60),ylim=c(0,5000))
dev.off()
jpeg("keep.proportion.jpeg")
hist(keep[,2],breaks=20,main="",xlab="Cumulation of hypermethylated samples",xlim=c(0,60),ylim=c(0,5000))
dev.off()



jpeg("keep.proportion.0.5.jpeg")
hist(keep[which(keep[,2]>34),2],breaks=20,main="",xlab="Cumulation of hypermethylated samples",col="red",lwd=1.5)
dev.off()

length(keep[which(keep[,2]>34),2])

gene248<-rownames(file1)[keep[which(keep[,2]>34),1]]
write.table(gene248,file="gene248.txt",row.names=F,col.names=F,quote=F)

gene248<-read.table("gene248.txt",sep="\t",as.is=T)

gene248anno<-read.table("gene248.anno.txt",sep="\t",as.is=T)
ncbigene<-read.table("gene_result.hsa.txt",head=T,sep="\t",as.is=T)

discover1<-ncbigene[na.omit(match(gene248anno[,7],ncbigene[,4])),]
write.table(discover1,file="discover1.txt",col.names=F,row.names=F,quote=F,sep="\t")


pdf("var.pdf")
hist(var,breaks=200)
dev.off()

d<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt.bed",sep="\t",as.is=T)

(0.063/0.309)-1
(30.733/182.9)-1





colnames(file1)







x<-rnorm(68,-2,1)
y<-rnorm(8,10,1)
wilcox.test(x,y)$p.value
qqplotR<-function(p){
  observed <- sort(p)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed)) 
  lexp <- -(log10(expected / (length(expected)+1)))
  data.frame(lobs,lexp)
}
pp1<-qqplotR(p)
pdf("rrbs.2.qqplot.pdf", width=6, height=6)
plot(c(0,8), c(0,8), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,8), ylim=c(0,8), las=1, xaxs="i", yaxs="i", bty="l")
points(pp1$lexp,pp1$lobs, pch=1, cex=.4, col="blue")   
abline(h=-log(0.05/nrow(pp1),10),lty=2,col="grey")
dev.off()

