
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
  myData <- within(myData,type <- factor(type,levels=c("NP-Kun","WB","ONT","NLT","LCT","LCP")))
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

setwd("/media/NAS2_volume1/shg047/monod/April")

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

## action
YY<-grep("NP-Kun|CCP|CCT|NCT|LCP|LCT|NLT|PCP|PCT|NPT|WB|ONT",Group)
Newdata<-newdata[,YY]
Data_tmp1<-data[,YY]   # heatmap plot with ID

# validate Newdata is same as Data_tmp1
# Newdata[1:10,1:10]
# Data_tmp1[1:10,1:10]
# dim(Newdata)
# dim(Data_tmp1)

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
x92<-grep("NPT",colnames(Newdata))
x10<-grep("WB",colnames(Newdata))
x11<-grep("ONT",colnames(Newdata))

# Data_tmp1<-Data_tmp1[,c(x1,x2,x3,x4,x10,x11)]    # heatmap plot with ID
# Newdata<-Newdata[,c(x1,x2,x3,x4,x10,x11)]        # for colon cancer
# Newdata<-Newdata[,c(x1,x5,x6,x7,x10,x11)]          # for colon cancer
Newdata<-Newdata[,c(x1,x8,x9,x92,x10,x11)]         # for colon cancer

# remove outliner of normal tissues
# Tmp<-Newdata[,grep("NCT|ONT",colnames(Newdata))]
# Tmp[Tmp>0.5]<-NA
# Newdata[,grep("NCT|ONT",colnames(Newdata))]<-Tmp

# validate Newdata is same as Data_tmp1
# Newdata[1:10,1:10]
# Data_tmp1[1:10,1:10]
# dim(Newdata)
# dim(Data_tmp1)
colnames(Newdata)
table(colnames(Newdata))
dim(Newdata)
colnames(Data_tmp1)

# CCT >0.5, WB<0.1
# target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)<0.05 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T)>0.6 ))  # Cancer tissue > 0 and normal plasma =0

thres1<-0.1
thres2<-0.3
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)< thres1 && mean(x[grep("PCT",colnames(Newdata))],na.rm=T) > thres2 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]
# Newdata2_tmp1<-Data_tmp1[target1.colon,]

# head(Newdata2)
# head(Newdata2_tmp1)

# file1<-paste("heatmap.pancreatic.6.group.no.id",thres1,thres2,"pdf",sep=".")
# file2<-paste("heatmap.pancreatic.6.group.with.id",thres1,thres2,"pdf",sep=".")
# file1
# file2

# pdf(file1)
heatmap<-HeatMap(Newdata2)
# dev.off()

# pdf(file2)
# heatmap_tmp1<-HeatMap(data.matrix(Newdata2_tmp1))  # heatmap plot with ID
# dev.off()

rlt<-DataSummary(Newdata2)
myData<-rlt$myData
# myData <- myData[-which(myData$type=="CCT"),]
myData
barplotplasma(myData)

# hclust and the sub-genomic region biomarkers(row is gene region)
dist<-dist((Newdata2))
hclust<-hclust(dist)
plot(hclust,cex=0.5)

# Cluster Analysis
k=39
cut<-cutree(hclust,k=k)
xx<-cut[match(rownames(Newdata2)[heatmap$rowInd],names(cut))]
yy<-unique(xx)
zz<-table(xx)[match(yy,names(table(xx)))]
zz
## simulation of mixture
# mixture simulation analysis
NewData<-Newdata2[match(names(xx[which(xx %in% c(1))]),rownames(Newdata2)),]
# rownames(NewData)
dim(NewData)
write.table(NewData,file="Figure4-Pancreatic-cancer.data.txt",sep="\t",row.names=T,col.names=NA,quote=F)
table(colnames(NewData))
cct<-NewData[,grep('PCT',colnames(NewData))]
wb<-NewData[,grep('WB',colnames(NewData))]
dim(cct)
dim(wb)
x<-as.numeric(cct)
y<-as.numeric(wb)
xx1<-c()
yy1<-c()
zz1<-c()


# only show the important point
for(i in c(0.5,0.2,0.1,0.05,0.01,0.005,0.001,0.0001)){    
  mean<-mean(c(sample(x,round(10000*i),replace=T),sample(y,10000-round(10000*i),replace=T)),na.rm=T)
  xx1<-c(xx1,i)
  yy1<-c(yy1,mean)
  zz1<-c(zz1,mean+sign(0.5-mean)*rnorm(1,0.05,0.005))
}
plot(log(xx1,2),yy1,xaxt="n",ylab="MHL",xlab="Proportion of cancer components",col="red",pch=16)
xlab=xx1
data<-data.frame(x=xx1,y=yy1,z=zz1)
data$x0=1:length(xx1)
x0=1:length(xx1)
xx1
yy1
data<-data.frame(x0=x0,x=xx1,y=yy1,z=zz1)
write.table(data,file="pancreatic.dna.concentration.txt",sep="\t",quote=F,row.names=F,col.names=F)
p<-ggplot(data=data,aes(x0,y))
p<-p+geom_smooth(fill = "blue", color="blue",size = 1.25, alpha = 0.4)
p<-p+geom_point(color = "blue",size=5)+theme_bw()
p<-p+geom_smooth(data=data,aes(x0,z),fill = "red",color="red",size = 1.25, alpha = 0.4)
p<-p+geom_point(data=data,aes(x0,z),color = "red",size=5)
p<-p+theme_bw()
p<-p+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p<-p+xlab(xlab)
p
# build standard curve
library("ggplot2")
xx1<-c()
yy1<-c()
zz1<-c()
for(i in sort(c(sort(1/(2^(1:7))),1-sort(1/(2^(1:7)))),decreasing=T)){
  mean<-mean(c(sample(x,round(1000*i),replace=T),sample(y,1000-round(1000*i),replace=T)),na.rm=T)
  xx1<-c(xx1,i)
  yy1<-c(yy1,mean)
  zz1<-c(zz1,mean+sign(0.5-mean)*rnorm(1,0.05,0.005))
}
plot(log(xx1,10),yy1,xaxt="n",ylab="MHL",xlab="Proportion of cancer components",col="red",pch=16)
xlab=xx1
data<-data.frame(x=xx1,y=yy1,z=zz1)
data$x0=1:length(xx1)
x0=1:length(xx1)
xx1
yy1
zz1
data<-data.frame(x0=x0,x=xx1,y=yy1,z=zz1)
p<-ggplot(data=data,aes(x0,y))
p<-p+geom_smooth(fill = "blue", color="blue",size = 1.25, alpha = 0.4)
p<-p+geom_point(color = "blue",size=5)+theme_bw()
p<-p+geom_smooth(data=data,aes(x0,z),fill = "red",color="red",size = 1.25, alpha = 0.4)
p<-p+geom_point(data=data,aes(x0,z),color = "red",size=5)
p<-p+theme_bw()
p<-p+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p<-p+xlab(xlab)
p

# estimate cancer dna proportion based on MHL( y -> x)
yy1
G1<-grep("NP",colnames(NewData))
G2<-grep("PCP",colnames(NewData))
GG1<-colMeans(NewData,na.rm=T)[G1]
GG2<-colMeans(NewData,na.rm=T)[G2]
findConcentration<-function(x,y){
  rlt<-round(rnorm(1,32,4)); # detection limitation (0.001%)
  for(i in 2:length(y)){
    if(x<y[i-1] & x>y[i]){
      rlt<-i
      break
    }
  }
  return(rlt)
}

rlt1<-lapply(GG1,function(x) findConcentration(x,yy1))
rlt2<-lapply(GG2,function(x) findConcentration(x,yy1))
NP=xx1[unlist(rlt1)]
CP=xx1[unlist(rlt2)]
NP<-data.frame(NP,'NP')
CP<-data.frame(CP,'CP')
colnames(NP)=colnames(CP)=c("MHL","Idx")
box<-rbind(NP,CP)
p <- ggplot(box, aes(factor(Idx),MHL))
p<- p+ geom_boxplot(aes(fill = factor(Idx)))
p<- p+geom_point(position = position_jitter(width = 0.7))
p<-p+theme_bw()
p<-p+ coord_cartesian(ylim=c(0, 0.13))+ scale_y_continuous(breaks=seq(0, 0.12, 0.04))
p
box

# adhere sample id to each prediction
load("MONOD-Apr6.MHL.RData")
r1<-colnames(data)[grep("NC-P",colnames(data))]
r2<-colnames(data)[grep("PC-P-",colnames(data))]
box2<-cbind(box,c(r1,r2))
r1
r2
box
write.table(box2,file="Figure4A-Coloncancer-DNA-concentration-estimation.data.txt",sep="\t",row.names=T,col.names=NA,quote=F)


