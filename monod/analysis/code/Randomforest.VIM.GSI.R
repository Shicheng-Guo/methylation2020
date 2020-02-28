rm(list=ls())
setwd("/home/sguo/Dropbox/Project/PanCancer")
load("rf.proximity.6631.1500tree.RData")
rlt3<-data.frame(rownames(RF$importance)[match(sort(RF$importance,decreasing=T)[1:300],RF$importance)],sort(RF$importance,decreasing=T)[1:300])
inf<-read.table("/home/sguo/monod/data/paad/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/jhu-usc.edu_PAAD.HumanMethylation450.9.lvl-3.TCGA-S4-A8RM-01A-11D-A378-05.txt",head=F,skip=2,sep="\t")
gene<-inf[match(rlt3[,1],inf[,1]),3]
chr<-inf[match(rlt3[,1],inf[,1]),4]
start<-inf[match(rlt3[,1],inf[,1]),5]-60
end<-inf[match(rlt3[,1],inf[,1]),5]+60
GSI<-gsi[match(rlt3[,1],names(gsi))]
cpg<-rlt3[,1]
VIM<-rlt3[,2]
rlt4<-data.frame(cpg,VIM,GSI,gene,chr,start,end)
write.table(rlt4,file="top300.txt",sep="\t",quote=F,col.names=NA,row.names=T)


