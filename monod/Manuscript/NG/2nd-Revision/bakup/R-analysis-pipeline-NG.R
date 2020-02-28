##===============================================================================================================================================================
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
## Don't provided barplot, Journal usually require boxplot
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
##   Can we just use one pipeline and run it and then get all the Figures? Why? (No, impossible, different analysis require different files and start points)
##   We can prepare code for each Figure, Tables and Supplementary Tables. 
##================================================================================================================================================================

## Figure 1

# collect all the functions, that's enough. 


## Figure 2

#!/usr/bin/R
# Methylation Haploinfo to Methylation haplotype load (mode weigth i)
# Contact: Shicheng Guo
# Version 1.3
# Update: May/29/2015

hap_comb_tmp<-function(hapinfo.single){
  # hap_com.single is a character of haplotype, such as "TCTCTT"
  hap_com.single<-c()
  for(i in 1:nchar(hapinfo.single)){
    for(j in i:(nchar(hapinfo.single))){
      hap_com.single<-c(hap_com.single,(substr(hapinfo.single,i,j)))
    } 
  }
  return(hap_com.single)
}

hap_comb<-function(hapinfo){
  # haplotype is array of observed methylation haplotype and return
  hap_comb<-unlist(lapply(hapinfo,hap_comb_tmp))
  return(hap_comb)
}

mhl<-function(hap_comb_array){
  mhl_value<-c()
  hap_comb_array<-unlist(hap_comb_array)
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  unmeth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="T"
  hap_comb_array<-c(hap_comb_array[meth_hap],hap_comb_array[unmeth_hap])
  for(i in unique(nchar(hap_comb_array))){
    input_tmp<-which(nchar(hap_comb_array)==i)  
    nmeth<-length(grep("C",hap_comb_array[input_tmp]))
    nunmeth<-length(grep("T",hap_comb_array[input_tmp]))
    mhl_value<-c(mhl_value,i*nmeth/(nunmeth+nmeth))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}

mf<-function(hapinfo){
  mf_value<-c()
  hap_comb_array<-unlist(hapinfo)
  meth_hap=unlist(lapply(hap_comb_array,function(x) unlist(strsplit(x,""))))
  nmeth<-length(grep("C",meth_hap))
  nunmeth<-length(grep("T",meth_hap))
  mf_value<-nmeth/(nunmeth+nmeth)
  return(mf_value)
}

# Figuure 2A
hapinfo<-c(rep("TTTT",16))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2B
hapinfo<-c(rep("CCCC",16))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2C
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
MHL<-mhl(hap_comb(unique(hapinfo)))
MF<-mf(hapinfo)
me<-methentropy(hapinfo)
# Figuure 2D
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2E
hapinfo<-c(rep("CCTT",8),rep("TTCC",8))
MHL<-mhl(hap_comb(unique(hapinfo)))

## Figure 3


## Figure 4
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

findConcentration<-function(x,y){
  rlt<-round(rnorm(1,35,2)); # detection limitation (0.001%)
  for(i in 2:length(y)){
    if(x<y[i-1] & x>y[i]){
      rlt<-i
      break
    }
  }
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
Newdata<-Newdata[,c(x1,x5,x6,x7,x10,x11)]          # for lung cancer
# Newdata<-Newdata[,c(x1,x8,x9,x92,x10,x11)]         # for pancreatic cancer

table(colnames(Newdata))
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

thres1<-0.05
thres2<-0.5
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)< thres1 && mean(x[grep("LCT",colnames(Newdata))],na.rm=T) > thres2 ))  # Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]

# head(Newdata2)
# head(Newdata2_tmp1)

file1<-paste("heatmap.lung.6.group.no.id",thres1,thres2,"pdf",sep=".")
file2<-paste("heatmap.lung.6.group.with.id",thres1,thres2,"pdf",sep=".")
file1
file2

pdf(file1)
heatmap<-HeatMap(Newdata2)
dev.off()

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

k=80
cut<-cutree(hclust,k=k)
cut
xx<-cut[match(rownames(Newdata2)[heatmap$rowInd],names(cut))]
xx
yy<-unique(xx)
yy
# sample size in each subgroup
zz<-table(xx)[match(yy,names(table(xx)))]
zz

# Group 2: 
Newdata2<-Newdata2[match(names(xx[which(xx %in% c(15))]),rownames(Newdata2)),]
input.long<-data.frame(variable=names(colMeans(Newdata2,na.rm = T)),value=colMeans(Newdata2,na.rm = T))
table(input.long$variable)
input.long <- within(input.long,variable <- factor(variable,levels=c("NP-Kun","WB","ONT","NLT","LCT","LCP")))
ggplot(aes(y = value, x = variable), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw() + ylim(0,1)+ theme(
  plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank()
)


## Figure 5

normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

#################################################
################## CCP ##########################
#################################################
par(mfrow=c(2,2))
hist(CCP[2,],breaks=30,xlim=c(0,60),col="darkblue",border="darkblue",xlab="Colon",ylab="Counts of cs-MHL",main="Colon")
hist(CCP[11,],breaks=10,xlim=c(0,60),col="darkblue",border="darkblue",xlab="CT",ylab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)

Zmax<-matrix(nrow=nrow(CCP),ncol=ncol(CCP))
for(i in 1:ncol(CCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,60)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    zmp <- (CCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(CCP)
colnames(Zmax)=colnames(CCP)
Zmax

par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i])
}

par(mfrow=c(2,2))
library("beeswarm")
input<-list(CCP=Zmax[2,],WB=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:4,pch=16,method="center",ylab="Z-Score")


p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="bonferroni"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)

sum(p1[2,]>-log(0.05,10))
length(unique(c(which(p1[2,]>-log(0.05,10)),which(p1[11,]>-log(0.05,10)))))

library("gplots")
heatmap.2(p1,Rowv=F,Colv=F,density.info="none",trace="none")


write.table(p1,file="colon.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)
getwd()

# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==2))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))/ncol(Zmax)

### Figure 
par(mfrow=c(4,3),mar=c(4,4,2,1))
AUC<-c()
for(k in 1:11){
  xx<-c()
  sen<-c()
  spe<-c()
  for(fdr in seq(0,1,by=0.01)){
    xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
    sen<-c(sen,xx/30)
    spe<-c(spe,1-fdr)
  }
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main=rownames(Zmax)[k])
  AUC<-c(AUC,mean(sen))
}

Zmax[k,]<-Zmax[2,]+Zmax[11,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Colon+CT")

# barplot of AUC for CRC cancer
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Colon+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")


#################################################
################## LCP ##########################
#################################################

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

Zmax<-matrix(nrow=nrow(LCP),ncol=ncol(LCP))
for(i in 1:ncol(LCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,29)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    
    zmp <- (LCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}

rownames(Zmax)=rownames(LCP)
colnames(Zmax)=colnames(LCP)
Zmax

par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i])
}

library("beeswarm")
input<-list(LCP=Zmax[6,],WB=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:4,pch=16,method="center",ylab="Z-Score")


par(mfrow=c(2,2))
hist(LCP[6,],breaks=20,xlim=c(0,60),col="darkblue",border="darkblue",ylab="Counts of samples",xlab="Counts of ls-MHL",main="Colon")
hist(LCP[11,],breaks=8,xlim=c(0,60),col="darkblue",border="darkblue",ylab="Counts of samples",xlab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)


# bar plot show Z-score for one lung sample

barplot(Zmax[,2],col="darkblue",ylab="Z-score")

file=list.file("*txt")
for(i in file){
  d<-read.table(i)
}

# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==6))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==6)))/ncol(Zmax)


p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="bonferroni"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)

length(unique(c(which(p1[6,]>-log(0.001,10)),which(p1[11,]>-log(0.001,10)))))

write.table(p1,file="lung.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)

## Figure C
par(mfrow=c(4,3),mar=c(4,4,2,1))
AUC<-c()
for(k in 1:11){
  xx<-c()
  sen<-c()
  spe<-c()
  for(fdr in seq(0,1,by=0.01)){
    xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
    sen<-c(sen,xx/30)
    spe<-c(spe,1-fdr)
  }
  AUC<-c(AUC,mean(sen))
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main=rownames(Zmax)[k])
}
Zmax[k,]<-Zmax[6,]+Zmax[11,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Lung+CT")

# barplot of AUC for Lung cancer
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Lung+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")


#################################################
################## PCP ##########################
#################################################
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))
PCP<-data.matrix(read.table("PCP.txt",sep="\t",head=T,row.names=1,as.is=T))

CCP<-PCP
par(mfrow=c(2,2))
hist(CCP[2,],breaks=30,xlim=c(0,60),col="darkblue",border="darkblue",xlab="Colon",ylab="Counts of cs-MHL",main="Colon")
hist(CCP[11,],breaks=10,xlim=c(0,60),col="darkblue",border="darkblue",xlab="CT",ylab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)

Zmax<-matrix(nrow=nrow(CCP),ncol=ncol(CCP))
for(i in 1:ncol(CCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,60)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    zmp <- (CCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(CCP)
colnames(Zmax)=colnames(CCP)
Zmax

par(mfrow=c(2,5),mar=c(4,2,4,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i],col="darkblue",ylim=c(-3,3))
}

par(mfrow=c(2,2))
library("beeswarm")
input<-list(CCP=Zmax[2,],WB=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:4,pch=16,method="center",ylab="Z-Score")


p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="bonferroni"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)

sum(p1[2,]>-log(0.05,10))
length(unique(c(which(p1[2,]>-log(0.05,10)),which(p1[11,]>-log(0.05,10)))))

library("gplots")
heatmap.2(p1,Rowv=F,Colv=F,density.info="none",trace="none")


write.table(p1,file="pancreatic.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)
getwd()

# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==2))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))/ncol(Zmax)

### Figure 
par(mfrow=c(4,3),mar=c(4,4,2,1))
AUC<-c()
for(k in 1:11){
  xx<-c()
  sen<-c()
  spe<-c()
  for(fdr in seq(0,1,by=0.01)){
    xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
    sen<-c(sen,xx/10)
    spe<-c(spe,1-fdr)
  }
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main=rownames(Zmax)[k])
  AUC<-c(AUC,mean(sen))
}

Zmax[k,]<-Zmax[2,]+Zmax[11,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/10)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Colon+CT")

# barplot of AUC for CRC cancer
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Colon+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")



####################################################
################## Normal ##########################
####################################################

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))
PCP<-data.matrix(read.table("PCP.txt",sep="\t",head=T,row.names=1,as.is=T))

par(mfrow=c(3,4))
for(i in 1:(nrow(normal))){
  hist(normal[i,],col="darkblue",breaks=30,xlim=c(0,50),xlab=rownames(normal)[i],border ="darkblue",main="",ylab="Counts of rsMHL",cex.axis=1.15,cex.lab=1.5)  
  abline(v=30,lty=3,col="red",lwd=2)
}
par(mfrow=c(3,4),mar=c(5,3,2,1))
for(i in 1:(nrow(normal))){
  hist(normal[i,],col="darkblue",breaks=30,xlim=c(0,50),xlab=rownames(normal)[i],border ="darkblue",main="",ylab="Counts of rsMHL",cex.axis=1.15,cex.lab=1.5)  
}

Zmax<-matrix(nrow=nrow(normal),ncol=ncol(normal))
for(i in 1:ncol(normal)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,75)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    
    zmp <- (normal[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}

rownames(Zmax)=rownames(normal)
colnames(Zmax)=colnames(normal)
Zmax



p<-pnorm(-abs(Zmax))
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1
p1<--log(p1,10)

sum(p[2,]<0.05)
sum(p[6,]<0.05)

write.table(p1,file="normal.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)


####################################################
################## CCP.Z-score.txt ##########################
####################################################


Zmax<-read.table("CCP.Z-score.txt")
par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i],col = "darkblue")
}

Zmax<-read.table("LCP.Z-score.txt")
par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i],col = "darkblue")
}

Zmax<-read.table("PCP.Z-score.txt")
par(mfrow=c(2,5),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i],col = "darkblue")
}

