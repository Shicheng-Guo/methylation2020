setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
load("MONOD-march30.MHL.RData")
sum(is.na(data))/(nrow(data)*ncol(data))
# re-group/collect the samples #colon cancer
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
table(colnames(newdata))

## For Colon Cancer
YY<-grep("NP-Kun|CCP",Group)
Newdata<-newdata[,YY]
dim(Newdata)
idx<-sapply(colnames(Newdata),function(x) unlist(strsplit(x,"[.]"))[1])
colnames(Newdata)<-idx
Newdata<-data.matrix(Newdata)
head(Newdata)
colnames(Newdata)
dim(Newdata)
Newdata<-Newdata[,c(grep("CCP",colnames(Newdata)),grep("NP-Kun",colnames(Newdata)))]
table(colnames(Newdata))
dim(Newdata)

NumNANP<-apply(Newdata[,31:56],1,function(x) sum(is.na(x)))
NumNACP<-apply(Newdata[,1:30],1,function(x) sum(is.na(x)))
par(mfrow=c(3,3))
hist(NumNANP,col="blue",cex.main=0.85,main=c("Normal Plasma"),xlab="Numbers of missing value")
hist(NumNACP,col="blue",cex.main=0.85,main=c("Cancer Plasma"),xlab="Numbers of missing value")
colnames(Newdata)
mydata<-Newdata
mydata<-RawNARemove(mydata)
dim(mydata)
mydata<-ColNARemove(mydata)
dim(mydata)

## For Lung Cancer
YY<-grep("NP-Kun|LCP",Group)
Newdata<-newdata[,YY]
dim(Newdata)
idx<-sapply(colnames(Newdata),function(x) unlist(strsplit(x,"[.]"))[1])
colnames(Newdata)<-idx
Newdata<-data.matrix(Newdata)
head(Newdata)
colnames(Newdata)
dim(Newdata)
Newdata<-Newdata[,c(grep("LCP",colnames(Newdata)),grep("NP-Kun",colnames(Newdata)))]
table(colnames(Newdata))
dim(Newdata)

NumNACP<-apply(Newdata[,1:29],1,function(x) sum(is.na(x)))
NumNANP<-apply(Newdata[,31:55],1,function(x) sum(is.na(x)))
hist(NumNANP,col="blue",cex.main=0.85,main=c("Normal Plasma"),xlab="Numbers of missing value")
hist(NumNACP,col="blue",cex.main=0.85,main=c("Cancer Plasma"),xlab="Numbers of missing value")
colnames(Newdata)
mydata<-Newdata
mydata<-RawNARemove(mydata)
dim(mydata)
mydata<-ColNARemove(mydata)
dim(mydata)
nrow(mydata)/nrow(newdata)

## For Pancreatic Cancer
YY<-grep("NP-Kun|PCP",Group)
Newdata<-newdata[,YY]
dim(Newdata)
idx<-sapply(colnames(Newdata),function(x) unlist(strsplit(x,"[.]"))[1])
colnames(Newdata)<-idx
Newdata<-data.matrix(Newdata)
head(Newdata)
colnames(Newdata)
dim(Newdata)
Newdata<-Newdata[,c(grep("PCP",colnames(Newdata)),grep("NP-Kun",colnames(Newdata)))]
table(colnames(Newdata))
dim(Newdata)

NumNACP<-apply(Newdata[,1:10],1,function(x) sum(is.na(x)))
NumNANP<-apply(Newdata[,11:26],1,function(x) sum(is.na(x)))
par(mfrow=c(3,3))
hist(NumNANP,col="blue",cex.main=0.85,main=c("Normal Plasma"),xlab="Numbers of missing value")
hist(NumNACP,col="blue",cex.main=0.85,main=c("Cancer Plasma"),xlab="Numbers of missing value")
colnames(Newdata)
mydata<-Newdata
mydata<-RawNARemove(mydata)
dim(mydata)
mydata<-ColNARemove(mydata)
dim(mydata)






