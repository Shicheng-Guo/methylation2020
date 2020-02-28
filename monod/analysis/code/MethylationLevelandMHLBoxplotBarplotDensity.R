---
title: "MF and MHL"
author: "Shicheng Guo"
date: "Monday, March 14, 2016"
output: html_document
---

Here is for boxplot and density plot for methylation level(the later will be for MHL). 

```{r}
data<-read.table("MergePancancerMethMatrix.txt",sep="\t",head=T,row.names=1,check.names=F)
data=data.matrix(data)
colnames(data)
xx<-strsplit(colnames(data)[grep("-T-|-P-",colnames(data))],"-")
y1<-unlist(lapply(xx,function(x) x[1]))
y2<-unlist(lapply(xx,function(x) x[2]))
y3<-unlist(lapply(xx,function(x) sprintf("%02d", as.numeric(x[3]))))

colnames(data)[grep("-T-|-P-",colnames(data))]<-paste(y1,y2,y3,sep="-")
newdata<-data.frame(value=as.numeric(data),Group=rep(colnames(data),each=nrow(data)))

myData <- aggregate(newdata$value,by =list(type=newdata$Group),
                    FUN = function(x) c(mean = mean(x,na.rm=T), 
                                        sd = sd(x,na.rm=T),
                                        sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                        me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
)

myData <- do.call(data.frame, myData)
colnames(myData)=c("type","mean","sd","sem","me")
# myData$type <- factor(myData$type, levels = c("NT","NP","TT","CP"))
myData
# myData$sd<-c(myData$sd[1],0.02,myData$sd[3],0.14)
# Plot one standard error (standard error of the mean/SEM)

library("ggplot2")
pdf("Value.barplot.ggplot2.pdf")
ylab="Sample Id"
xlab="Methylation Level"
title="Methylation Level Distribution"
ggplot(myData, aes(x =type, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue",width=0.8) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),size=1.0) +
  ggtitle(title) + 
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  coord_flip()+
  xlab(xlab) +
  ylim(0,1)+
  ylab(ylab)+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=5),axis.text.y = element_text(hjust=0))
dev.off()

pdf("Value.boxplot.ggplot2.pdf")
ggplot(newdata, aes(factor(Group),value)) + 
  geom_boxplot(aes(fill = factor(Group)),outlier.shape=NA)+ 
  coord_flip()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=5),axis.text.y = element_text(hjust=0),legend.position="none")
dev.off()
ggsave("Value.boxplot.ggplot2.pdf")


pdf("CRC-P-01.pdf")
par(mfrow=c(2,1))
ncol=match("CRC-P-01",colnames(data))
plot(density(data[,ncol],na.rm=T),lty=1,col=1)
ncol=match("CRC-P-11",colnames(data))
lines(density(data[,ncol],na.rm=T),lty=2,col=2)
legend("topright",legend=c("CRC-P-01","CRC-P-11"),lty=c(1,2),col=c(1,2))
dev.off()
```

Here is for boxplot and density plot for MHL level. 

```{r}
data<-read.table("/oasis/tscc/scratch/shg047/monod/mhl/mhl-03-18-2016.txt",sep="\t",head=T,row.names=1,check.names=F)
data=data.matrix(data)
colnames(data)
xx<-strsplit(colnames(data)[grep("-T-|-P-",colnames(data))],"-")
y1<-unlist(lapply(xx,function(x) x[1]))
y2<-unlist(lapply(xx,function(x) x[2]))
y3<-unlist(lapply(xx,function(x) sprintf("%02d", as.numeric(x[3]))))

colnames(data)[grep("-T-|-P-",colnames(data))]<-paste(y1,y2,y3,sep="-")
newdata<-data.frame(value=as.numeric(data),Group=rep(colnames(data),each=nrow(data)))

myData <- aggregate(newdata$value,by =list(type=newdata$Group),
                    FUN = function(x) c(mean = mean(x,na.rm=T), 
                                        sd = sd(x,na.rm=T),
                                        sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                        me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
)

myData <- do.call(data.frame, myData)
colnames(myData)=c("type","mean","sd","sem","me")
# myData$type <- factor(myData$type, levels = c("NT","NP","TT","CP"))
myData
# myData$sd<-c(myData$sd[1],0.02,myData$sd[3],0.14)
# Plot one standard error (standard error of the mean/SEM)

library("ggplot2")
pdf("MHL.barplot.ggplot2.pdf")
ylab="Sample Id"
xlab="Methylation Level"
title="Methylation Level Distribution"
ggplot(myData, aes(x =type, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue",width=0.8) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),size=1.0) +
  ggtitle(title) + 
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  coord_flip()+
  xlab(xlab) +
  ylim(0,1)+
  ylab(ylab)+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=5),axis.text.y = element_text(hjust=0))
dev.off()

pdf("MHL.boxplot.ggplot2.pdf")
ggplot(newdata, aes(factor(Group),value)) + 
  geom_boxplot(aes(fill = factor(Group)),outlier.shape=NA)+ 
  coord_flip()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=5),axis.text.y = element_text(hjust=0),legend.position="none")
dev.off()
ggsave("Value.boxplot.ggplot2.pdf")

pdf("6-P-01.MHL.pdf")
par(mfrow=c(2,1))
ncol=match("6-P-11",colnames(data))
plot(density(data[,ncol],na.rm=T),lty=2,col=2)
ncol=match("6-P-01",colnames(data))
lines(density(data[,ncol],na.rm=T),lty=1,col=1)
legend("topright",legend=c("CRC-P-01","CRC-P-11"),lty=c(1,2),col=c(1,2))
dev.off()

```
