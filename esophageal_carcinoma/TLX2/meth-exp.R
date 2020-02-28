setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/esca/tlx2")
data<-read.table("meth-mRNA.txt",row.names=1,head=T)

summary(lm(data$me5~log(data$mRNA,10)))
pdf("ZNF.meth.exp.pdf")
par(mfrow=c(2,2))
plot(log(data$EXP) ~ MF, data = data,type="n",col=as.numeric(as.factor(unlist(lapply(as.character(data$Target),function(x) substr(x,1,1)))))+1,xlab="Average Methylation Level")
points(y=log(data$EXP),x=data$MF,pch=as.numeric(as.factor(unlist(lapply(as.character(data$Target),function(x) substr(x,1,1)))))+1)
abline(lm(log(data$EXP)~data$MF),lwd=2,lty=5) # regression line (y~x) 
legend("topright",legend=c("Normal","Cancer"),pch=c(2,3),bty="n")
dev.off()

