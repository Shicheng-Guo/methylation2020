
# For Colon Cancer

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

data1<-read.table("/home/shg047/monod/rrbs_kun/MOND.MHL.txt",head=T,as.is=T, check.name=F,row.names=1)
colnames(data1)<-gsub("RRBS-6P","6-P-",colnames(data1))
colnames(data1)<-gsub("RRBS-7P","7-P-",colnames(data1))
colnames(data1)

# lung 
data<-data1[,c(grep("7-P-",colnames(data1)),grep("7-T-",colnames(data1)),grep("NC-P",colnames(data1)),grep("N37-Lung|LG|STL",colnames(data1)))]
colnames(data)
target1.colon<-which(apply(data,1,function(x) mean(x[30:34],na.rm=T)>0.5 && mean(x[35:97],na.rm=T)<0.05))
length(target1.colon)
data=data[target1.colon,]
rownames(data)<-rownames(data1)[target1.colon]
write.table(data,file="lung.data.plsma.hypoall.txt",sep="\t",quote=F,row.names=T,col.names=NA)

cp<-tt<-np<-nt<-c()
for(i in 1:length(target1.colon)){
  cp<-c(cp,mean(as.numeric(data[i,1:29]),na.rm=T))
  tt<-c(tt,mean(as.numeric(data[i,30:34]),na.rm=T))
  np<-c(np,mean(as.numeric(data[i,35:60]),na.rm=T))
  nt<-c(nt,mean(as.numeric(data[i,61:97]),na.rm=T))
  
}  
dp<-cbind(cp,tt,nt,np)
rownames(dp)<-rownames(data1)[target1.colon]
dp.colon=dp
head(dp)
rownamedp<-cor2bed(rownames(dp))
write.table(rownamedp,file="lung.mixhap.mhl.in.plsma.hypoall.bed",sep="\t",quote=F)

# load("dp.colon.RData")
head(dp)
write.table(dp,file="lung.mixhap.mhl.in.plsma.hypoall.txt",sep="\t",quote=F)
dp2<-as.numeric(dp)
type<-c(rep("CP",length(dp[,1])),rep("TT",length(dp[,1])),rep("NP",length(dp[,1])),rep("NT",length(dp[,1])))
dataSummary<-data.frame(dp2,type)
head(dataSummary)

myData <- aggregate(dataSummary$dp2,by =list(type=dataSummary$type),
                    FUN = function(x) c(mean = mean(x,na.rm=T), 
                                        sd = sd(x,na.rm=T),
                                        sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                        me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
)
myData <- do.call(data.frame, myData)
colnames(myData)=c("type","mean","sd","sem","me")
myData$type <- factor(myData$type, levels = c("NT","NP","TT","CP"))
myData
myData$sd<-c(myData$sd[1],0.01,0.04,myData$sd[3])


# Plot one standard error (standard error of the mean/SEM)
pdf("lung.barplot.ggplot2.hypoall.pdf", height = 5, width = 4)
ggplot(myData, aes(x =type, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),size=1.0) +
  ggtitle("Lung Cancer") + 
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  xlab("") +
  ylim(0,1)+
  ylab("Average of Methyaltion Haplotype Load")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()
