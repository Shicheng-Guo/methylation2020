
DataSummary<-function(data){
  rlt<-list()
  NewData<-data.frame(value=as.numeric(data),Group=rep(colnames(data),each=nrow(data)))
  ttest<-t.test(NewData[grep("CP",NewData$Group),1],NewData[grep("Kun",NewData$Group),1],na.rm=T)
  myData <- aggregate(NewData$value,by =list(type=NewData$Group),
                      FUN = function(x) c(mean = mean(x,na.rm=T), 
                                          sd = sd(x,na.rm=T),
                                          sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                          me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
  )
  myData <- do.call(data.frame, myData)
  colnames(myData)=c("type","mean","sd","sem","me")
  rlt$ttest<-ttest
  rlt$myData<-myData
  return(rlt)
}

HeatMap<-function(data){
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  ColSideColors=sidecol(colnames(data))
  heatmap.2(data,trace="none",cexRow = 0.1,cexCol = 0.7, ColSideColors=ColSideColors,density.info="none",col=colors,Colv=F,keysize=0.9, margins = c(5, 10))
}

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")

load("MONOD-Apr6.MHL.RData")
## Miss Value Statistic
sum(is.na(data))/(nrow(data)*ncol(data))
na<-apply(data,2,function(x) sum(is.na(x)))
sort(na)
barplot(sort(na),horiz=T,las=1,cex.names=0.25,col="blue")
write.table(sort(na),file="missing.value.table.txt",sep="\t")

## re-group/collect the samples #colon cancer
Group<-colnames(data)
# Group[grep("CTR|NC-|Pregn",Group)]<-"NP"
Group[grep("NC-",Group)]<-"NP-Kun"
Group[grep("CTR|Pregn",Group)]<-"NP-Dennis"
Group[grep("6-P-",Group)]<-"CCP"
Group[grep("6-T-|HCT116|Colon_Tumor_Primary|CTT-|metastasis_colon|tumor_colon",Group)]<-"CCT"
Group[grep("normal_colon|N37-Colon|SG-01",Group)]<-"NCT"

Group[grep("7-P-",Group)]<-"LCP"
Group[grep("7-T-|adenocarcinoma_lung|tumor_lung",Group)]<-"LCT"
Group[grep("LG-01|N37-Lung|normal_lung",Group)]<-"NLT"

Group[grep("PC-P-",Group)]<-"PCP"
Group[grep("PC-T-",Group)]<-"PCT"
Group[grep("PA-01|N37-Pancreas",Group)]<-"NPT"

Group[grep("WB-",Group)]<-"WB"
Group[grep("STL|N37-|methylC-",Group)]<-"ONT"
newdata<-data
colnames(newdata)<-Group
newdata<-data.frame(newdata,check.names=F)

## action
YY<-grep("NP-Kun|CCP|CCT|NCT|LCP|LCT|NLT|PCP|PCT|NPT|WB|ONT",Group)
Newdata<-newdata[,YY]

idx<-sapply(colnames(Newdata),function(x) unlist(strsplit(x,"[.]"))[1])
colnames(Newdata)<-idx
Newdata<-data.matrix(Newdata)
head(Newdata)
colnames(Newdata)
table(colnames(Newdata))
dim(Newdata)

x1<-grep("NP-Kun",colnames(Newdata))
x2<-grep("CCP",colnames(Newdata))
x3<-grep("CCT",colnames(Newdata))
x4<-grep("NCT",colnames(Newdata))
x5<-grep("LCP",colnames(Newdata))
x6<-grep("LCT",colnames(Newdata))
x7<-grep("NLT",colnames(Newdata))
x8<-grep("PCP",colnames(Newdata))
x9<-grep("PCT",colnames(Newdata))
x10<-grep("WB",colnames(Newdata))
x11<-grep("ONT",colnames(Newdata))

Newdata<-Newdata[,c(x1,x2,x3,x4,x10,x11)]
colnames(Newdata)
table(colnames(Newdata))
dim(Newdata)

# only CCP high, all the value were high
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)<0.2 && mean(x[grep("CCP",colnames(Newdata))],na.rm=T)>0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.5))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
DataSummary(Newdata2)
heatmap(Newdata2)

# only CCP high, all the value were high
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)<0.2 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.3))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
DataSummary(Newdata2)
heatmap(Newdata2)

# only CCT high, all the value were high
target1.colon<-which(apply(Newdata,1,function(x)   mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.3 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
DataSummary(Newdata2)
# heatmap(Newdata2)

# CCT high while NCT, ONT is low
target1.colon<-which(apply(Newdata,1,function(x) mean(x[grep("NCT",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.3 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
rlt<-DataSummary(Newdata2)
# heatmap<-HeatMap(Newdata2)
myData<-rlt$myData
myData
# pdf("xxxx.pdf")
barplotplasma(myData)
# dev.off()

# CCT >0.5, WB<0.1
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.5 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
rlt<-DataSummary(Newdata2)
heatmap<-HeatMap(Newdata2)
myData<-rlt$myData
myData

# WB<0.1, NCT<0.1, CCT>0.5, ONT<0.1
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("NCT",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.5 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)<0.1))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
rlt<-DataSummary(Newdata2)
rlt
myData<-rlt$myData
myData
heatmap<-HeatMap(Newdata2)
dim(Newdata2)

dist<-dist((Newdata2))
hclust<-hclust(dist)
plot(hclust,cex=0.5)

# CCT high while WB, NP low 
target1.colon<-which(apply(Newdata,1,function(x) mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("NCT",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.3))
length(target1.colon)
Newdata2=Newdata[target1.colon,]
rlt<-DataSummary(Newdata2)
rlt
myData<-rlt$myData
myData
barplotplasma(myData)

heatmap<-HeatMap(Newdata2)

dist<-dist((Newdata2))
hclust<-hclust(dist)
plot(hclust,cex=0.5)

cut<-cutree(hclust,k=5)
xx<-cut[match(rownames(Newdata2)[heatmap$rowInd],names(cut))]
yy<-unique(xx)
yy

length(which(xx==yy[1]))
length(which(xx==yy[2]))
length(which(xx==yy[3]))
length(which(xx==yy[4]))
length(which(xx==yy[5]))


NewData<-Newdata2[match(names(xx[which(xx==yy[1])]),rownames(Newdata2)),]
rlt<-DataSummary(NewData)
rlt

NewData<-Newdata2[match(names(xx[which(xx==yy[2])]),rownames(Newdata2)),]
rlt<-DataSummary(NewData)
myData<-rlt$myData
pdf("group-2.pdf")
barplotplasma(myData)
dev.off()

NewData<-Newdata2[match(names(xx[which(xx==yy[3])]),rownames(Newdata2)),]
rlt<-DataSummary(NewData)
myData<-rlt$myData
myData
pdf("group-3.pdf")
barplotplasma(myData)
dev.off()

NewData<-Newdata2[match(names(xx[which(xx==yy[4]|xx==yy[5])]),rownames(Newdata2)),]
rlt<-DataSummary(NewData)
myData<-rlt$myData
pdf("group-5.pdf")
barplotplasma(myData)
dev.off()

# CCT high while WB, NP low 
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("NP",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.7 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
DataSummary(Newdata2)
HeatMap(Newdata2)

# CCT + while WB, ONT low 
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
DataSummary(Newdata2)
HeatMap(Newdata2)

# CCT + while WB, ONT low 
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("NCT",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)<0.1 && any(x[grep("CCT",colnames(Newdata))]>0.3)))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
rlt<-DataSummary(Newdata2)
rlt
# heatmap(Newdata2)


# CCT high while WB, NP low 
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.1 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)<0.8 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.3 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
rlt<-DataSummary(Newdata2)
rlt
# heatmap(Newdata2)

#
target1.colon<-which(apply(Newdata,1,function(x) mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.3 && mean(x[grep("HNT",colnames(Newdata))],na.rm=T)<0.2 && mean(x[grep("NP",colnames(Newdata))],na.rm=T)))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)

# write.table(NewData,file="Supplementary-colon.cancer.tissue.hyper.txt",col.names=NA,row.names=T,quote=F,sep="\t")

myData<-rlt$myData

barplotplasma<-function(myData){
library("ggplot2")
# pdf("Lung.cancer.mhl.plasma.groups.107.pdf")
ylab="Average MHL"
xlab="Group"
title="Average MHL in different Groups"
myData
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
