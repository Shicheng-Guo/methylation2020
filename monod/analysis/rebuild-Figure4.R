
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
  
#  note: this function include correlation based heatmap (pearson or spearman)
#  data: row is gene and column is sample
#  colname and rowname cannot be NULL  
#  Usage example:
#  test<- matrix(runif(100),nrow=20)
#  colnames(test)=c("A","A","A","B","B")
#  rownames(test)=paste("Gene",1:20,sep="")
#  HeatMap(test)
  
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  
  Hclust=function(x){hclust(x,method="complete")}
  
  Distfun=function(x){as.dist((1 - cor(t(x),method = "spearman")))}
  
  ColSideColors=sidecol(colnames(data))
  
  heatmap.2(data,trace="none", 
            hclust=Hclust,
            distfun=Distfun, 
            cexRow = 1, cexCol = 1,
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=T,Rowv = TRUE,
            keysize=0.9, margins = c(5, 10)
            )
}



setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
load("MONOD-Apr6.MHL.RData")
## Miss Value Statistic
# sum(is.na(data))/(nrow(data)*ncol(data))
# na<-apply(data,2,function(x) sum(is.na(x)))
# sort(na)
# barplot(sort(na),horiz=T,las=1,cex.names=0.25,col="blue")
#write.table(sort(na),file="missing.value.table.txt",sep="\t")

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
x10<-grep("WB",colnames(Newdata))
x11<-grep("ONT",colnames(Newdata))

Data_tmp1<-Data_tmp1[,c(x1,x2,x3,x4,x10,x11)]  # heatmap plot with ID
Newdata<-Newdata[,c(x1,x2,x3,x4,x10,x11)]

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
thres2<-0.5
target1.colon<-which(apply(Newdata,1,function(x) mean(x[grep("NCT",colnames(Newdata))],na.rm=T)< 0.1 && mean(x[grep("ONT",colnames(Newdata))],na.rm=T)< 0.1 && mean(x[grep("WB",colnames(Newdata))],na.rm=T)< 0.1 && mean(x[grep("CCT",colnames(Newdata))],na.rm=T) > thres2))  # Cancer tissue > 0 and normal plasma =0

length(target1.colon)
Newdata2=Newdata[target1.colon,]
head(Newdata2)

input.long<-data.frame(variable=names(colMeans(Newdata2,na.rm = T)),value=colMeans(Newdata2,na.rm = T))
input.long <- within(input.long,variable <- factor(variable,levels=c("NP-Kun","WB","ONT","NCT","CCT","CCP")))
ggplot(aes(y = value, x = variable), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw() + ylim(0,1)+ theme(
  plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank()
) 

rlt<-DataSummary(Newdata2)
myData<-rlt$myData
barplotplasma(myData)

barplotplasma<-function(myData){
  library("ggplot2")
  # pdf("Lung.cancer.mhl.plasma.groups.107.pdf")
  ylab="Average MHL"
  xlab="Group"
  title="Average MHL in different Groups"
  myData <- within(myData,type <- factor(type,levels=c("NP-Kun","WB","ONT","NCT","CCT","CCP")))
  library("ggplot2")
  ggplot(myData, aes(x =type, y = mean)) +  
    geom_bar(position = position_dodge(), stat="identity",width=0.9) + 
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.0) +
    ggtitle(title) + 
    theme_bw() +
    theme(panel.grid.major = element_blank())+
    xlab(xlab) +
    ylab(ylab)+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10),axis.text.y = element_text(hjust=0))+
    geom_jitter()
}

library("reshape2")

Row<-c()
for(i in 1:100){
NewData2<-Newdata2
NewData2<-NewData2[sample(1:nrow(NewData2),500,replace = T),]
id=1:nrow(NewData2)
input<-data.frame(id,NewData2)
colnames(input)<-c("id",unlist(lapply(colnames(NewData2),function(x) unlist(strsplit(x,"[.]"))[1])))
input.long<-melt(input, id.vars=c("id"))
Row<-rbind(Row,tapply(input.long$value,input.long$variable,function(x) mean(x,na.rm=T)))
}

input<-data.frame(id=1:nrow(Row),Row)
input.long<-melt(input,id.vars=c("id"))
head(input.long)
table(input.long$variable)
input.long<-within(input.long,variable<-factor(variable,levels = c("NP.Kun","WB","ONT","NCT","CCT","CCP")))
ggplot(aes(y = value, x = variable), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw() +  theme(
  plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank()
) 

barplotplasma(myData)

