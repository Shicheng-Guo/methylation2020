
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\HMP")
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis")
data=read.table("HMP.txt",head=T,sep="\t")
data=read.table("mhl.ecf.txt",head=T,sep="\t")
set.seed(111)
head(data)
data=na.omit(data)
data[,3]<-data[,3]+abs(rnorm(nrow(data),0.001,0.01))
data<-subset(data,data[,2]<100 & data[,3]<0.1 )
head(data)
data[,1]
plot(data[,3],data[,2],col=as.numeric(as.factor(data[,1]))+1)
fit<-lm(data[,2]~data[,3])
fit
summary(fit)
abline(fit,col="red")
plot(y=data$Yield,x=data$MHL,pch=19,cex=1,col=as.numeric(as.factor(data$Lab.ID))+1,ylab="Normalized yield from 1ml plasma",xlab="Estimated tumor fraction")
legend("topright",legend=c("CP","LP","NP"),col=c(2,3,4),pch=19)
fit<-lm(Yield~MHL,data)
abline(fit,lwd=3,lty=2,col="brown")

summary(fit)
write.table(data,file="last.update.adjust.txt",quote=F,row.names=T)


Lung cancer  29
Colon cancer	30
Normal plasma	75

x<-c(29,30,75)
names(x)<-c("Lung cancer","Colon cancer","Normal plasma")
pie(x,col=rainbow(3),clockwise=TRUE)


pdf("pie.pdf")
par(mfrow=c(2,2))
x<-c(0.652122693,0.120601817,0.059117375,0.168158116)
names(x)<-c("WB","CT","NT","TBA")
pie(x,col=c("red","green","blue","NA"),clockwise=TRUE,main="Cancer Plasma",labels=round(x,2))

x<-c(0.706402936,0.069151898,0.22)
names(x)<-c("WB","NT","TBA")
pie(x,col=c("red","blue","NA"),clockwise=TRUE,main="Normal Plasma",labels=round(x,2))
dev.off()


pdf("pie-2.pdf")
library(plotrix)

par(mfrow=c(2,2))
x<-c(0.652122693,0.120601817,0.059117375,0.168158116)
names(x)<-c("WB","CT","NT","TBA")
pie(x,col=c("red","green","blue","NA"),clockwise=TRUE,main="Cancer Plasma",labels=round(x,2))

x<-c(0.706402936,0.069151898,0.22)
names(x)<-c("WB","NT","TBA")
pie(x,col=c("red","blue","NA"),clockwise=TRUE,main="Normal Plasma",labels=round(x,2))
dev.off()



# Normal Plasma from Dennis
data<-read.table("Dennis.Normal2.MF",head=F,sep="\t",row.names=1)
rowMean<-rowMeans(data,na.rm = T)
hypo<-which(rowMean<0.3)
pdf("AML.in.MHB.Dennis.pdf")
hist(rowMean)
dev.off()

# 








