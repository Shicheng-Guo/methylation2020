

setwd("/home/sguo/monod/rrbs")
pos<-c()
file<-list.files(pattern="*.trim")
for (i in 1:length(file)){
  tmp<-read.table(file[i],sep="\t",as.is=T)
  pos<-c(pos,tmp[,1])
  print(i)
}


encode.rrbs.pos.freq<-sort(table(pos))
save(encode.rrbs.pos.freq,file="encode.rrbs.pos.freq.RData")



setwd("G:\\monod\\encode")
load("encode.rrbs.pos.freq.RData")
ls()
head(encode.rrbs.pos.freq)
pdf("encode.rrbs.pos.freq.2.pdf",width=2.5,height=2.5)
hist(encode.rrbs.pos.freq,breaks=150,main="",xlab="")
dev.off()


setwd("G:\\monod\\encode")
load("RRBS.MethylationBlock.RData")

pdf("RRBS.encode.methyblock.pearson.correlation.distribution.pdf",width=4.05,height=4.36)
hist(rrbsmb$hdrc[,1],breaks=200,col="red",border="red",main="",xlab="Pearson correlation",cex.lab=1.2,cex.axis=1.2,lwd=0.1)
dev.off()

hist(rrbsmb$hdr[,3],breaks=100,xlim=c(1,50))
par(mar=c(4,2,1,1))
hist(rrbsmb$hdr[,3],breaks=150,xlim=c(1,44),ylim=c(0,8000),col="red",border="red",main="",xlab="",ylab="",cex.lab=1.1,cex.axis=0.70,lwd=0.1,xaxt="n",yaxt="n")
axis(1,at=seq(4,55,by=10),labels=seq(4,55,by=10),cex.axis=1.2)
axis(2,at=c(0,2,4,6,8,10,12)*1000,labels=c(0,2,4,6,8,10,12),pos=4,cex.axis=1.2)
title(xlab="Pearson correlation",line =2.5,cex.lab=1.25)
title(ylab = "Counts(10K)",line = 0.5,cex.lab=1.25)

pdf("Figure.S2.length.of.hdr.pdf")
hist(rrbsmb$hdr[,3],breaks=150,xlim=c(1,44),ylim=c(0,8000),col="red",border="red",main="",xlab="",ylab="",cex.lab=1.1,cex.axis=0.70,lwd=0.1,yaxt="n",xaxt="n")
axis(1,at=seq(4,55,by=10),labels=seq(4,55,by=10),cex.axis=1.2)
axis(2,at=c(0,2,4,6,8,10,12)*1000,labels=c(0,2,4,6,8,10,12),pos=4,cex.axis=1.2)
title(xlab="Length of methylation blocks",line =2.5,cex.lab=1.25)
title(ylab = "Counts(10K)",line = 0.5,cex.lab=1.25)
dev.off()

