#!/usr/bin/R
setwd("/home/gsc/Dropbox/Project/methylation/array/ld")
install.packages("skmeans")

load("sk2.RData")
load("sk3.RData")
load("sk4.RData")
load("sk5.RData")
load("sk6.RData")
load("sk7.RData")
load("sk8.RData")
load("sk9.RData")

Prediction<-cbind(sparty$cluster,sparty3$cluster,sk4$cluster,sk5$cluster,sk6$cluster,sk7$cluster,sk8$cluster,sk9$cluster)
load("OvClinic.RData")
ClinicNumericFact<-matrix(NA,412,8)
for (i in 2:9){
  ClinicNumericFact[,i-1]<-as.numeric(clinc[,i])
}
colnames(ClinicNumericFact)<-tolower(colnames(clinc)[2:9])

data1<-cbind(ClinicNumericFact,Prediction)

if(length(which(apply(data1,1,function(x) sum(is.na(x)))>0))>0){
  data2<-data1[-which(apply(data1,1,function(x) sum(is.na(x)))>0),]
}else{
  data2<-data1;
}
dim(data2)

dim(data1)[1]-dim(data2)[1]
ClinicNumericFact<-data2[,1:8]
Prediction<-data2[,9:16]

CorrectRand<-Vi<-matrix(NA,8,8)

for (i in 1:dim(ClinicNumericFact)[2]){
  for (j in 1:dim(Prediction)[2]){
     result<-cluster.stats(NULL,ClinicNumericFact[,i],Prediction[,j],compareonly=T)
     CorrectRand[i,j]<-result$corrected.rand
     Vi[i,j]<-result$vi
  }
}

color2D.matplot(Vi,show.values=T,show.legend=T,xlab="Cluster By FPCA",ylab="Clinical Cluster",main="2D matrix plot")
color2D.matplot(CorrectRand,show.values=3,show.legend=T,xlab="Cluster By FPCA",ylab="Clinical Cluster",main="CorrectRand 2D matrix plot")

va<-cbind(clinc$TUMORGRADE,sparty3$cluster)
xtabs(~va[,1]+va[,2])

# the clust effect is not very good, However, what does array work?
# merge+cluser+compare

setwd("/home/sguo/score")
load("ChrChrosomeName.cor.RData")
block<-array()
j<-0
for (i in seq(11,(ncol(cor)-10),by=5)){
  j<-j+1
  block[j]<-mean(cor[(i-10):i,i:(i+10)],na.rm=T)
}

  which(block>0.75)
  save(block,file="ChrChrosomeName.block.RData")







