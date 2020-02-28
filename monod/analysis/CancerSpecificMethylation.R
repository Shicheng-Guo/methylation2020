
# group specificity index,GSI
setwd("../pan")
load("PanPairMethData.RData")
pheno=data$pheno
xmean <- rowsum(data[,2:ncol(data)], pheno)/table(pheno)
gsi<-apply(xmean,2,function(x) (length(x)-sum(x)/max(x))/(length(x)-1))
rlt<-data.frame(pheno=names(table(pheno)),xmean[,match(names(sort(gsi,decreasing=T)[1:5]),colnames(xmean))])
write.table(rlt,file="pancancer.high.gsi.site.txt",sep="\t",quote=F,col.names=NA,row.names=T)

jpeg("gsi.distribution.jpeg")
hist(gsi,xlab="Histogram of Group Specific Index",main="")
dev.off()

num90<-sum(gsi>0.85)
xmean90<-xmean[,match(names(sort(gsi,decreasing=T)[1:num90]),colnames(xmean))]
xx1<-table(data$pheno)[as.numeric(names(table(apply(xmean90,2,function(x) which.max(x)))))]
xx2<-table(apply(xmean90,2,function(x) which.max(x)))
rlt2<-data.frame(sample_size=xx1,xx2)
write.table(rlt2,file="target.status.specifi.txt",sep="\t",quote=F,col.names=NA,row.names=T)

gsi2<-gsi[order(gsi,decreasing=T)]
rlt3<-c()
for(i in seq(1,22,by=2)){
  z<-0
  j<-0
  while(j<length(gsi2) & z<5){
  j=j+1
  max<-which.max(xmean[,match(names(gsi2[j]),colnames(xmean))])
  if(max==i){
    z<-z+1
    rlt3<-rbind(rlt3,(c(names(gsi2[j]),names(table(pheno))[i])))
  }
}
}

rlt3<-data.frame(rlt3)
inf<-read.table("/home/sguo/monod/data/paad/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/jhu-usc.edu_PAAD.HumanMethylation450.9.lvl-3.TCGA-S4-A8RM-01A-11D-A378-05.txt",head=F,skip=2,sep="\t")
cpg<-rlt3[,1]
cancer<-substr(rlt3[,2],1,4)
gene<-inf[match(rlt3[,1],inf[,1]),3]
chr<-inf[match(rlt3[,1],inf[,1]),4]
start<-inf[match(rlt3[,1],inf[,1]),5]-60
end<-inf[match(rlt3[,1],inf[,1]),5]+60
GSI<-gsi[match(rlt3[,1],names(gsi))]
rlt4<-data.frame(cpg,cancer,GSI,gene,chr,start,end)
write.table(rlt4,file="target.status.specifi.site.txt",sep="\t",quote=F,col.names=NA,row.names=T)



