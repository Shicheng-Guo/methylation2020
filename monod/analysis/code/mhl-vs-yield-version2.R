setwd("/home/sguo/Dropbox/Project/methylation/monod/analysis/code")
data<-read.table("mhl-vs-yiled.txt")  # Supplementary Table
plot(data[,8],data[,4])
data<-subset(data,data[,4]<0.12 & data[,4]>0 & data[,8]<120 & data[,8]>0) # local alignment
plot(data[,4],data[,8])
summary(lm(data[,8]~data[,4]))
plot(log(data[,8],10),log(data[,4],10))
summary(lm(log(data[,4],10)~log(data[,8],10)))

pdf("mhl-vs-yield.pdf",width=5,height=5)
plot(data[,4],data[,8],type="n",ylim=c(0,110),xlab="Estimated tumor DNA fraction",xlim=c(0,0.1),ylab="Normalized yield from 1ml plasma")
npc<-which(data[,3]=="NP")
lc<-which(data[,3]=="LCP")
cc<-which(data[,3]=="CCP")
points(y=data[npc,8],x=data[npc,4],pch=16,cex=2,col="blue")
points(y=data[lc,8],x=data[lc,4],pch=16,cex=2,col="green")
points(y=data[cc,8],x=data[cc,4],pch=16,cex=2,col="red")
abline(lm(data[,8]~data[,4]),lty=5,lwd=4,col="brown")
dev.off()
cancer<-subset(data,data[,3]=="CCP" |data[,3]=="LCP" )
normal<-subset(data,data[,3]=="NP")
fit1<-lm(cancer[,8]~cancer[,4])
fit1
fit2<-lm(normal[,8]~normal[,4])
fit2
