#! Remember all the advance analysis should be conducted in Genome-miner so that everything can be documented.
# system("scp /oasis/tscc/scratch/shg047/monod/hapinfo/N37Salk.saminfo shg047@genome-miner.ucsd.edu:/media/NAS1/shg047/monod/hapinfo")
# system("scp /media/NAS1/shg047/tcga/pancancer/colon.450h.RData shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/monod/hapinfo")

library("monod")
library("ggplot2")
library("reshape2")

setwd("/media/NAS1/shg047/monod/hapinfo")
system("scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/saminfo/N37Salk.saminfo ./")
system("scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/hapinfo/MHL4.RData ./")

load("/media/NAS1/shg047/monod/hapinfo/MHL4.RData")
data<-read.table("/media/Home_Raid1/shg047/NAS1/mhl/mhl.v1.4.merged",head=T,row.names=1)
saminfo<-read.table("/media/NAS1/shg047/monod/hapinfo/N37Salk.saminfo",sep="\t")

Data=data[,grep("STL|N37|age|ENC|new|centenarian|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
Data<-rename(Data)
colnames(Data)<-unlist(lapply(colnames(Data),function(x) unlist(strsplit(x,"[.|-]"))[1]))
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
table(colnames(DATA))
gsirlt<-gsi(DATA)
save(gsirlt,file="gsirlt.RData")

# validate the biomarkers
library("monod")
load("/media/NAS1/shg047/tcga/pancancer/PancancerMethMatrix_March2016.RData")
load("../hapinfo/gsirlt.RData")
colon<-subset(gsirlt,group=="Colon")

# make sure colon specfici markers from WGBS are hypemethyation 
cg<-bed2cg(cor2bed(colon[,1]))
input<-data[unique(match(cg[,8],rownames(data))),grep("COAD",colnames(data))]
input<-input[which(apply(input,1,function(x) mean(x,na.rm=T))>0.2),]
bed<-unique(rbedintersect(cor2bed(colon[,1]),unique(cg2bed(rownames(input))[,1:4]))[,1:4])

# validation performance
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC|CCP|LCP|NP",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
test<-DATA[match(bed[,4],rownames(DATA)),]
wbcexclude1<-which(apply(test[,grep("WBC",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude2<-which(apply(test[,grep("Spleen",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude3<-which(apply(test[,grep("Stomach",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude4<-which(apply(test[,grep("Kidney",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude5<-which(apply(test[,grep("Intestine",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude6<-which(apply(test[,grep("Liver",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude7<-which(apply(test[,grep("Pancreas",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude8<-which(apply(test[,grep("Lung",colnames(test))],1,function(x) mean(x,na.rm=T))>0.0001)
wbcexclude9<-which(apply(test[,grep("Brain",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude<-unique(c(wbcexclude1,wbcexclude2,wbcexclude3,wbcexclude4,wbcexclude5,wbcexclude6,wbcexclude7,wbcexclude8,wbcexclude9))
wbcexclude<-unique(c(wbcexclude1))
test<-test[-wbcexclude,]
xx<-apply(test,2,function(x) sum(na.omit(x)>0.1)/length(na.omit(x)))
xx<-data.frame(mdf=xx,type=names(xx))
xx$type <- factor(xx$type,levels = c("NP","LCP","CCP","WBC","Brain","Colon","Intestine",
                                     "Kidney","Liver","Lung","Pancreas","Spleen","Stomach"
),ordered = TRUE)
uxx<-melt(test)
colnames(uxx)<-c("type","mhl")
uxx$type <- factor(uxx$type,levels = c("NP","LCP","CCP","WBC","Brain","Colon","Intestine",
                                       "Kidney","Liver","Lung","Pancreas","Spleen","Stomach"
),ordered = TRUE)
png("assign.number.png")
p1<-ggplot(aes(y = mhl, x = type, fill = type), data = uxx) + 
  geom_boxplot(outlier.shape =NA,outlier.colour="blue")+ 
  coord_flip()+
  geom_point(position = position_jitter(width = 0.2),size=0.5)

p2<-ggplot(aes(y = mdf, x = type, fill = type), data = xx) + 
  geom_boxplot(outlier.shape =NA,outlier.colour="blue")+ 
  coord_flip()+
  geom_point(position = position_jitter(width = 0.2),size=0.5)
multiplot(p1, p2, cols=2)
dev.off()
