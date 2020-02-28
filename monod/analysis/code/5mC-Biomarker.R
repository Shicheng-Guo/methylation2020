


##############################################################################################
########## Biomarker identification based on 5mC of MONOD dataset ############################
##############################################################################################
setwd("/home/shg047/oasis/monod/BAM/MF_PileOMeth")
data=read.table("project.rlt.txt",head=T,row.names=1,sep="\t",as.is=T,check.names=F)
CRCP<-data[,grep("CRC-P",colnames(data))]
CRCT<-data[,grep("CRC-T",colnames(data))]
LCP<-data[,grep("LC-P",colnames(data))]
LCT<-data[,grep("LC-T",colnames(data))]
PCP<-data[,grep("PC-P",colnames(data))]
PCT<-data[,grep("PC-T",colnames(data))]
NCP<-data[,grep("NC-P",colnames(data))]
dim(CRCP)
dim(CRCT)
dim(LCP)
dim(LCT)
dim(PCP)
dim(PCT)
dim(NCP)
idx<-unlist(lapply(strsplit(colnames(data),"-0"),function(x) x[[1]]))
newdata=data
colnames(newdata)=idx
rlt<-data.frame(t(apply(newdata,1,function(x) tapply(x,idx,function(x) mean(x,na.rm=T)))))
rlt<-data.frame(t(rlt))
newrlt1<-subset(rlt,NC.P<0.2 & LC.T>LC.P & LC.P>0)
newrlt2<-subset(rlt,NC.P<0.2 & CRC.T>CRC.P & CRC.P>0)
newrlt3<-subset(rlt,NC.P<0.2 & PC.T>PC.P & PC.P>0)
result<-unique(rbind(newrlt1,newrlt2,newrlt3))
for(i in rownames(result)){
len=length(na.omit(as.numeric(data[grep(i,rownames(data)),])))
print(len)
}
write.table(result,file="Biomarker-based-on-5mC-CRC-LC-NCP.txt",col.names=NA,row.names=T,sep="\t",quote=F)


##############################################################################################
########## Biomarker identification based on MHL of MONOD dataset ############################
##############################################################################################
setwd("/home/shg047/oasis/monod/hapinfo/June")
saminfo<-read.table("/home/shg047/oasis/monod/saminfo.txt",sep="\t")
# data<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/June/monod.mhl.june25.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
# save(data,file="monod.mhl.june25.RData")
load("monod.mhl.june25.RData")
colnames(data)[grep("STL",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("STL",colnames(data))],saminfo[,1]),2])
colnames(data)[grep("WB",colnames(data))]<-"WBC"
colnames(data)[grep("N37",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("N37",colnames(data))],saminfo[,1]),2])
colnames(data)[grep("methylC",colnames(data))]<-"H1"
colnames(data)[grep("6-T|Colon_Tumor_Primary",colnames(data))]<-"CCT"
colnames(data)[grep("7-T-",colnames(data))]<-"LCT"
colnames(data)[grep("6-P-",colnames(data))]<-"CCP"
colnames(data)[grep("7-P-",colnames(data))]<-"LCP"
colnames(data)[grep("NC-P-",colnames(data))]<-"NCP"
colnames(data)[grep("PC-P-",colnames(data))]<-"PCP"
colnames(data)[grep("PC-T-",colnames(data))]<-"PCT"
data=data[,-grep("mix",colnames(data))]
data<-data[,grep("CCP|CCT|LCP|LCT|PCP|PCT|NCP",colnames(data))]
idx<-unlist(lapply(strsplit(colnames(data),"[.]"),function(x) x[[1]]))

newdata=data
colnames(newdata)=idx
rlt<-data.frame(t(apply(newdata,1,function(x) tapply(x,idx,function(x) mean(x,na.rm=T)))))
newrlt1<-subset(rlt,NCP<0.01 & LCT>LCP & LCP>0.1)
newrlt2<-subset(rlt,NCP<0.01 & CCT>CCP & CCP>0.1)
newrlt3<-subset(rlt,NCP<0.01 & PCT>PCP & PCP>0.1)
result<-unique(rbind(newrlt1,newrlt2,newrlt3))
dim(result)
write.table(result,file="Biomarker-based-on-MHL-CRC-LC-PC-NCP.txt",col.names=NA,row.names=T,sep="\t",quote=F)

########################################################################################################
############################## Overlap Biomarker Regions ###############################################
########################################################################################################
r1<-read.table("/home/shg047/oasis/monod/BAM/MF_PileOMeth/Biomarker-based-on-5mC-CRC-LC-NCP.txt",head=T,row.names=1,sep="\t",as.is=T)
r2<-read.table("/home/shg047/oasis/monod/hapinfo/June/Biomarker-based-on-MHL-CRC-LC-PC-NCP.txt",head=T,row.names=1,sep="\t",as.is=T)
r1[which(rownames(r1) %in% rownames(r2)),]
r2[which(rownames(r2) %in% rownames(r1)),]











