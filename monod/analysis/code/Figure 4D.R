---
  title: "cancer speicfic haplotype and mhl in tumor plasma.Rmd"
author: "Shicheng Guo"
date: "Monday, February 01, 2016"
output: html_document
---
  Update Figure 4C. 


pre-load function is as the following:
  ```{r}
RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
} 

gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,mean)))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(data[i,][which(index==group[j])]))/10^(mean(data[i,][which(index==gmax)])))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    gsi<-c(gsi,gsit)
    print(c(gmax,gsit))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi)
  return(rlt)
}

```


```{r}
data1<-read.table("/home/shg047/monod/rrbs_kun/rrbs.kun.mhl.na.impute.txt",head=T,as.is=T, check.name=F)
colnames(data1)<-gsub("RRBS-6P","6-P-",colnames(data1))
colnames(data1)<-gsub("RRBS-7P","7-P-",colnames(data1))
colnames(data1)

# colon 
data<-data1[,c(grep("6-P",colnames(data1)),grep("6-T",colnames(data1)),grep("NC-P",colnames(data1))),grep("[N37-Colon|STL001SG]",colnames(data))]
colnames(data)
target1.colon<-which(apply(data,1,function(x) sum(x[31:35]>0)>1 && sum(x[36:58]>0.2)==0))
length(target1.colon)
data=data[target1.colon,]
colon.data=data
cp<-tp<-np<-c()
for(i in 1:length(target1.colon)){
  cp<-c(cp,mean(as.numeric(data[i,1:30])))
  tp<-c(tp,mean(as.numeric(data[i,31:35])))
  np<-c(np,mean(as.numeric(data[i,36:55])))
}  
dp<-cbind(cp,tp,np)
rownames(dp)<-rownames(data)[target1.colon]
dp.colon=dp
# load("dp.colon.RData")
head(dp)
dp2<-as.numeric(dp)
type<-c(rep("CP",length(dp[,1])),rep("TP",length(dp[,1])),rep("NP",length(dp[,1])))
dataSummary<-data.frame(dp2,type)
head(dataSummary)

myData <- aggregate(dataSummary$dp2,by =list(type=dataSummary$type),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        sem=sd(x)/sqrt(length(x)),
                                        me=qt(1-0.05/2,df=length(x)*sd(x)/sqrt(length(x))))
)
myData <- do.call(data.frame, myData)
myData

colnames(myData)=c("type","mean","sd","n","sem","me")
myData
myData$type <- factor(myData$type, levels = c("CP","TP","NP"))
myData

# Plot one standard error (standard error of the mean/SEM)
pdf("barplot.ggplot2.pdf", height = 3, width = 3)
tiff("barplot.ggplot2.tiff",res=300)
ggplot(myData, aes(x =type, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.1) +
  ggtitle("Colon Cancer") + 
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  xlab("") +
  ylab("Average of Methyaltion Haplotype Load")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()

```
# lung cancer 
```{r}
data<-data1[,c(grep("7-P",colnames(data1)),grep("7-T",colnames(data1)),grep("NC-P",colnames(data1)))]
colnames(data)
target1.lung<-which(apply(data,1,function(x) sum(x[30:34]>0)==5 && sum(x[35:54]>0.01)==0))
length(target1.lung)
data=data[target1.lung,]
cp<-tp<-np<-c()
for(i in 1:length(target1)){
  cp<-c(cp,mean(as.numeric(data[i,1:29])))
  tp<-c(tp,mean(as.numeric(data[i,30:34])))
  np<-c(np,mean(as.numeric(data[i,35:54])))
}  
dp<-cbind(cp,tp,np)
rownames(dp)<-rownames(data)[target1.lung]
dp.lung=dp

save(dp,file="dp.lung.RData")
load("dp.lung.RData")
head(dp)
dp2<-as.numeric(dp)
type<-c(rep("CP",length(dp[,1])),rep("TP",length(dp[,1])),rep("NP",length(dp[,1])))
dataSummary<-data.frame(dp2,type)
head(dataSummary)

myData <- aggregate(dataSummary$dp2,by =list(type=dataSummary$type),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        n =as.numeric(length(x)),
                                        sem=sd(x)/sqrt(n),
                                        me=qt(1-0.05/2,df=length(x)*sd(x)/sqrt(n)))
)
myData <- do.call(data.frame, myData)
colnames(myData)=c("type","mean","sd","n","sem","me")
myData
myData$type <- factor(myData$type, levels = c("CP","TP","NP"))
myData

# Plot one standard error (standard error of the mean/SEM)
pdf("barplot.ggplot2.pdf", height = 3, width = 3)

tiff("lung.cancer.tiff",res=300)
ggplot(myData, aes(x =type, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.1) +
  ggtitle("Lung Cancer") + 
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  xlab("") +
  ylab("Average of Methyaltion Haplotype Load")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()

```
For lung cancer, 740 regions


```{r}
data<-data1[,c(grep("PC-P",colnames(data1)),grep("PC-T",colnames(data1)),grep("NC-P",colnames(data1)))]
colnames(data)
target1.pancreatic<-which(apply(data,1,function(x) sum(x[11:15]>0)==5 && sum(x[16:35]>0.01)==0))
length(target1.pancreatic)
data=data[target1.pancreatic,]
cp<-tp<-np<-c()
for(i in 1:length(target1)){
  cp<-c(cp,mean(as.numeric(data[i,1:10])))
  tp<-c(tp,mean(as.numeric(data[i,11:15])))
  np<-c(np,mean(as.numeric(data[i,16:35])))
}  
dp<-cbind(cp,tp,np)
rownames(dp)<-rownames(data)[target1.pancreatic]
dp.lung=dp


save(dp,file="dp.pancreatic.RData")
load("dp.pancreatic.RData")
head(dp)
dp2<-as.numeric(dp)
type<-c(rep("CP",length(dp[,1])),rep("TP",length(dp[,1])),rep("NP",length(dp[,1])))
dataSummary<-data.frame(dp2,type)
head(dataSummary)

myData <- aggregate(dataSummary$dp2,by =list(type=dataSummary$type),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        n =as.numeric(length(x)),
                                        sem=sd(x)/sqrt(n),
                                        me=qt(1-0.05/2,df=length(x)*sd(x)/sqrt(n)))
)
myData <- do.call(data.frame, myData)
colnames(myData)=c("type","mean","sd","n","sem","me")
myData
myData$type <- factor(myData$type, levels = c("CP","TP","NP"))
myData

# Plot one standard error (standard error of the mean/SEM)
pdf("barplot.ggplot2.pdf", height = 3, width = 3)

tiff("lung.cancer.tiff",res=300)
ggplot(myData, aes(x =type, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.15) +
  ggtitle("Lung Cancer") + 
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  xlab("") +
  ylab("Average of Methyaltion Haplotype Load")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()
```
