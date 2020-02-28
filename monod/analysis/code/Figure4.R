#!/usr/bin/R

# source("https://bioconductor.org/biocLite.R")
# install.packages("preprocessCore")
# library("preprocessCore")
# library(gplots)
# library("preprocessCore")

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

barplotplasma<-function(myData){
  library("ggplot2")
  # pdf("Lung.cancer.mhl.plasma.groups.107.pdf")
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


load("C:\\Users\\shg047\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016\\MONOD-Apr6.MHL.RData")
load("/media/Home_Raid1/shg047/work/monod/May/MONOD-Apr6.MHL.RData")
## Miss Value Statistic
# sum(is.na(data))/(nrow(data)*ncol(data))
# na<-apply(data,2,function(x) sum(is.na(x)))
# sort(na)
# barplot(sort(na),horiz=T,las=1,cex.names=0.25,col="blue")
#write.table(sort(na),file="missing.value.table.txt",sep="\t")

## re-group/collect the samples #colon cancer
Group<-colnames(data)
newdata<-data[,grep("STL|_h1|WB|N37|SRX381716|SRX381719|SRX381722|SRX381569|Colon_Tumor|SRX381713|SRX381553",Group)]
colnames(newdata)
Group<-colnames(newdata)

rename<-function(colname){
  colname[grep("6-T-|HCT116|Colon_Tumor_Primary|CTT-|metastasis_colon|tumor_colon",colname)]<-"CCT"
  colname[grep("normal_colon|N37-Colon|SG-01",colname)]<-"NCT"
  colname[grep("7-T-|adenocarcinoma_lung|tumor_lung",colname)]<-"LCT"
  colname[grep("LG-01|N37-Lung|normal_lung",colname)]<-"NLT"
  
  
}
newdata<-data.matrix(newdata[-which(apply(newdata,1,function(x) sum(is.na(x))/ncol(newdata))>0.3),])
head(newdata)
newdata2<-normalize.quantiles(newdata,copy=TRUE)
colnames(newdata)<-colnames(newdata)
input<-data.matrix(newdata[order(apply(newdata,1,function(x) var(x,na.rm=T)),decreasing=T)[1:1500],])
head(input)

# define most variable regions as tissue specific regions.

heatmap.2(input,
           col = colorpanel(100,"blue","green","yellow"),
           margins = c(12, 22),
           trace = "none", 
           xlab = "Comparison",
           lhei = c(2, 8),
           scale = c("none"),
           symbreaks = min(input, na.rm=TRUE),
           na.color="blue",
           cexRow = 0.5, cexCol = 0.7,
           main = "MHL genes", 
           dendrogram = "both", 
           Colv = T )


x1# Group[grep("CTR|NC-|Pregn",Group)]<-"NP"
Group[grep("NC-",Group)]<-"NP-Kun"
Group[grep("CTR|Pregn",Group)]<-"NP-Dennis"
Group[grep("6-P-",Group)]<-"CCP"
Group[grep("6-T-|HCT116|Colon_Tumor_Primary|CTT-|metastasis_colon|tumor_colon",Group)]<-"CCT"
Group[grep("normal_colon|N37-Colon|SG-01",Group)]<-"NCT"
# Lung cancer
Group[grep("7-P-",Group)]<-"LCP"
Group[grep("7-T-|adenocarcinoma_lung|tumor_lung",Group)]<-"LCT"
Group[grep("LG-01|N37-Lung|normal_lung",Group)]<-"NLT"
# Pancreatic cancer
Group[grep("PC-P-",Group)]<-"PCP"
Group[grep("PC-T-",Group)]<-"PCT"
Group[grep("PA-01|N37-Pancreas",Group)]<-"NPT"
# WB and Normal tissues
Group[grep("WB-",Group)]<-"WB"
Group[grep("STL|N37-|methylC-",Group)]<-"ONT"

# data proxy
newdata<-data
colnames(newdata)<-Group
newdata<-data.frame(newdata,check.names=F)

colnames(newdata)

## action
YY<-grep("NP-Kun|CCP|CCT|NCT|LCP|LCT|NLT|PCP|PCT|NPT|WB|ONT",Group)
Newdata<-newdata[,YY]
Data_tmp1<-data[,YY]   # heatmap plot with ID

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

Data_tmp1<-Data_tmp1[,c(x1,x2,x3,x4,x10,x11)]  # heatmap plot with ID
Newdata<-Newdata[,c(x1,x2,x3,x4,x10,x11)]

colnames(Newdata)
table(colnames(Newdata))
dim(Newdata)
colnames(Data_tmp1)

thres1<-0.1
thres2<-0.5
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)< thres1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T) > thres2 ))  # Cancer tissue > 0 and normal plasma =0

length(target1.colon)
Newdata2=Newdata[target1.colon,]
Newdata2_tmp1<-Data_tmp1[target1.colon,]
colnames(Newdata2_tmp1)

pdf(file1)
heatmap<-HeatMap(Newdata2)
dev.off()

rlt<-DataSummary(Newdata2)
myData<-rlt$myData

# hclust and the sub-genomic region biomarkers(row is gene region)
dist<-dist((Newdata2))
hclust<-hclust(dist)
plot(hclust,cex=0.5)

k=5
cut<-cutree(hclust,k=k)
cut
xx<-cut[match(rownames(Newdata2)[heatmap$rowInd],names(cut))]
yy<-unique(xx)
yy
# sample size in each subgroup
zz<-table(xx)[match(yy,names(table(xx)))]
zz
