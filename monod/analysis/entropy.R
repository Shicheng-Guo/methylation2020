
entropy<-c()
for(p in 1:10){
N=16
b=4
data=matrix(sample(c(0,1),b*N,replace=T),N,b)
ent<-c()
for(i in 1:N){
  tmp<-c()
  for(j in 1:b){
  tmp=paste(tmp,data[i,j],sep="")
  }
  ent<-c(ent,tmp)
}
entropy<-c(entropy,entropy(table(ent),unit="log2")/b)
}


setwd("C:/Users/shicheng/Dropbox/Project/methylation/monod/phase3")


mbs<-read.table("WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed")
head(mbs)

name1<-c("Intergenic","Downstream","Enhancer","Exon","Intron","Promoter","UTR3","UTR5","miRNA")
q1<-c(15446,3452,13797, 7579,41966,8358,2463,9000,92)
name2<-c("Downstream","Enhancer","Exon","Intron","Promoter","UTR3","UTR5","miRNA")
q2<-c(1084,2876,3329,3486,3232,646,1680,24)
names(q1)=name1
names(q2)=name2
  
q3<-c(15446,3452,13797, 7579,41966,8358,2463,9000,92)/sum(q1)
q4<-c(1084,2876,3329,3486,3232,646,1680,24)/sum(q2)
names(q3)=name1
names(q4)=name2
  
pdf("distribution1a.pdf")
par(mfrow=c(2,2),mar=c(1,2,3,2))
pie(q1,col=rainbow(length(q1)),lwd=2)
pie(q2,col=rainbow(length(q2)),lwd=2)
par(mar=c(5,6,1,3))
bpt1<-barplot(q1,col=rainbow(length(q3)),lwd=2,horiz=T,las=2)
text(x= q1+3000, y= bpt1, labels=as.character(q1), xpd=TRUE)
bpt2<-barplot(q2,col=rainbow(length(q4)),lwd=2,horiz=T,las=2)
text(x= q2+200, y= bpt2, labels=as.character(q2), xpd=TRUE)
dev.off()


pdf("distribution1b.pdf")
par(mfrow=c(2,2),mar=c(1,2,3,2))
pie(q1,col=rainbow(length(q1)),lwd=2)
pie(q2,col=rainbow(length(q2)),lwd=2)
par(mar=c(5,6,1,3))
bpt1<-barplot(q3,col=rainbow(length(q3)),lwd=2,horiz=T,las=2)
text(x= q3+0.03, y= bpt1, labels=round(q3,3), xpd=TRUE)
bpt2<-barplot(q4,col=rainbow(length(q4)),lwd=2,horiz=T,las=2)
text(x= q4+0.015, y= bpt2, labels=round(q4,3), xpd=TRUE)
dev.off()
