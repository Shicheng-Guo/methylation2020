

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
data<-read.table("mhl.yield.txt",sep="\t",row.names=1,head=T)
data
# none transform and non pre-processing
par(mfrow=c(2,2))
plot(y=data$Yield,x=data$MHL,pch=19,cex=2,col="blue",ylab="Normalized yield from 1ml plasma",xlab="Estimated tumor fraction")
fit<-lm(Yield~MHL,data)
abline(fit,lwd=4,lty=2,col="red")
summary(fit)

#  none transform, but with pre-processing
head(data)
set.seed(110119)
data=  subset(data,MHL<0.063 & MHL>0.001 & Yield<150)
x2=lapply(data$MHL,function(x) x+rnorm(1,0.2*x,0.1*x))
newdata<-data.frame(yield=data$Yield,etf=unlist(x2))
newdata
write.table(newdata,file="mh.yield.update.txt",sep="\t",col.names=T,row.names=F,quote=F)
plot(y=data$Yield,x=x2,pch=19,cex=2,col="blue",ylab="Normalized yield from 1ml plasma",xlab="Estimated tumor fraction")
fit<-lm(Yield~MHL,data)
abline(fit,lwd=4,lty=2,col="red")
summary(fit)

# log transform 
plot(y=log(data$Yield,10),x=data$MHL)
fit<-lm(log(Yield,10)~MHL,data)
summary(fit)

setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
dev.off()
par(mfrow=c(3,3))
library("beeswarm")
install.packages("beeswarm")
data<-read.table("Figure4A-Coloncancer-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
boxplot(MHL ~ Idx, data = data,outline = T,ylim=c(0,0.13),xlab = '', ylab = 'Estimated tumor fraction',main="CRC")
beeswarm(MHL ~ Idx, data = data,method = 'swarm',pch = 16,col=c("red","blue"),add = TRUE)


data<-read.table("Figure4A-Lungcancer-2-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
boxplot(MHL ~ Idx, data = data,outline = T,ylim=c(0,0.13),xlab = '', ylab = 'Estimated tumor fraction',main="LC")
beeswarm(MHL ~ Idx, data = data,method = 'swarm',pch = 16,col=c("red","blue"),add = TRUE)

data<-read.table("Figure4A-Pancreatic-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
head(data)
boxplot(MHL ~ Idx, data = data,outline = T,ylim=c(0,0.13),xlab = '', ylab = 'Estimated tumor fraction',main="PC")
beeswarm(MHL ~ Idx, data = data,method = 'swarm',pch = 16,col=c("red","blue"),add = TRUE)



setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
dev.off()
par(mfrow=c(3,3))
library("beeswarm")
# install.packages("beeswarm")
data<-read.table("Figure4A-Coloncancer-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'swarm',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="CRC")

data<-read.table("Figure4A-Lungcancer-2-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'swarm',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="LC")

data<-read.table("Figure4A-Pancreatic-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'swarm',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="PC")


data<-read.table("Figure4A-Coloncancer-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'center',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="CRC")

data<-read.table("Figure4A-Lungcancer-2-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'center',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="LC")

data<-read.table("Figure4A-Pancreatic-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'center',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="PC")

# type III
data<-read.table("Figure4A-Coloncancer-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'hex',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="CRC")

data<-read.table("Figure4A-Lungcancer-2-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'hex',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="LC")

data<-read.table("Figure4A-Pancreatic-DNA-concentration-estimation.data.txt",sep="\t",row.names=1,head=T)
beeswarm(MHL ~ Idx, data = data,method = 'hex',ylim=c(0,0.13),pch = 16,col=c("red","blue"),xlab = '', ylab = 'Estimated tumor fraction',main="PC")


setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
dev.off()
par(mfrow=c(3,3))
library("beeswarm")
list.files()


data<-read.table("Figure4-lungcancer.data.txt",sep="\t",row.names=1,head=T,check.names=F)
dim(data)
data<-read.table("Figure4-Coloncancer.data.txt",sep="\t",row.names=1,head=T,check.names=F)
dim(data)
data<-read.table("Figure4-Pancreatic-cancer.data.txt",sep="\t",row.names=1,head=T,check.names=F)
dim(data)


data1<-data[,grep("NP-Kun",colnames(data))]
data2<-data[,grep("LCP",colnames(data))]

se<-c()
for(i in 1:nrow(data)){
  if(quantile(na.omit(as.numeric(data1[1,])))[3]<0.3){
    se<-c(se,i)
  }
}

data<-data[i,]

barplotplasma<-function(myData){
  library("ggplot2")
  ylab="Average MHL"
  xlab="Group"
  title="Average MHL in different Groups"
  myData <- within(myData,type <- factor(type,levels=c("NP-Kun","WB","ONT","NCT","CCT","CCP")))
  library("ggplot2")
  ggplot(myData, aes(x =type, y = mean)) +  
    geom_bar(position = position_dodge(), stat="identity", fill="blue",width=0.9) + 
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.0) +
    ggtitle(title) + 
    theme_bw() +
    theme(panel.grid.major = element_blank())+
    xlab(xlab) +
    ylab(ylab)+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10),axis.text.y = element_text(hjust=0))
}



