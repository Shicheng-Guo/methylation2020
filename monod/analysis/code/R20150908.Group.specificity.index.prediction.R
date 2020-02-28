
data<-read.table(file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
# heat only have 3 high GSI regions

file1<-read.table("RRBS_methHap_load_matrix_July2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file1)
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
remove=c("6-T-3","6-T-4","7-T-2",paste("NC-P-",19:24,sep=""),"PC-P-10","6-P-6",paste("PC-P-",c(2,3,6,9),sep=""))
file1<-file1[,-match(remove,samplename2)]
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
new<-read.table("saminfo.txt",sep="\t",as.is=T)
cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
matrix=file1
samplename2<-gsub("6-P","CC-P",samplename2)
samplename2<-gsub("7-P","LC-P",samplename2)
samplename2<-gsub("6-T","CC-T",samplename2)
samplename2<-gsub("7-T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
colnames(matrix)=samplename2

lung.signature<-subset(data,GSI>0.689 & group=="Lung")   # number=10
colon.signature<-subset(data,GSI>0.7709 & group=="Colon") # number=10
pancrease.signature<-subset(data,GSI>0.8474 & group=="Pancreas") # number=10

lung.signature<-subset(data,GSI>0.66 & group=="Lung")
colon.signature<-subset(data,GSI>0.75 & group=="Colon")
pancrease.signature<-subset(data,GSI>0.79 & group=="Pancreas")

lung.signature<-subset(data,GSI>0.66 & group=="Lung")
colon.signature<-subset(data,GSI>0.75 & group=="Colon")
pancrease.signature<-subset(data,GSI>0.83 & group=="Pancreas")

nrow(lung.signature) 
nrow(colon.signature)
nrow(pancrease.signature)

bed1<-cor2bed(rownames(matrix))
bed2<-cor2bed(lung.signature[,1])
bed3<-cor2bed(colon.signature[,1])
bed4<-cor2bed(pancrease.signature[,1])

rlt.lung<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed2,opt.string="-wa -u")
rlt.colon<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed3,opt.string="-wa -u")
rlt.pancrease<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed4,opt.string="-wa -u")

cor.lung<-bed2cor(rlt.lung)
cor.colon<-bed2cor(rlt.colon)
cor.pancrease<-bed2cor(rlt.pancrease)

length(cor.lung)
length(cor.colon)
length(cor.pancrease)

data.prediction.lung<-file1[match(cor.lung,rownames(file1)),c(grep("6-P",colnames(file1)))]
data.prediction.colon<-file1[match(cor.lung,rownames(file1)),c(grep("7-P",colnames(file1)))]
data.prediction.pancreatic<-file1[match(cor.lung,rownames(file1)),c(grep("PC-P",colnames(file1)))]

t1<-apply(data.prediction.lung,2,function(x) sum(x>0)/length(x))
t2<-apply(data.prediction.colon,2,function(x) sum(x>0)/length(x))
t3<-apply(data.prediction.pancreatic,2,function(x) sum(x>0)/length(x))


#
choose.prediction<-unique(c(cor.lung,cor.colon,cor.pancrease))
data.prediction.data<-file1[match(choose.prediction,rownames(file1)),c(grep("6-P",colnames(file1)),grep("7-P",colnames(file1)),grep("PC-P",colnames(file1)))]
t<-apply(data.prediction.data,2,function(x) { u=x>0;print(c(sum(unlist(u)[1:10]),sum(unlist(u)[11:20]),sum(unlist(u)[21:30])))})





