install.packages("~/monod_1.1.tar.gz")
library("monod")
library("ggplot2")
library("reshape2")

setwd("/media/NAS1/shg047/monod/hapinfo")
#system("scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/saminfo/N37Salk.saminfo ./")
#system("scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/hapinfo/MHL4.RData ./")
#data1<-read.table("/media/Home_Raid1/shg047/NAS1/monod/mhl/mhl.mhbs",head=T,row.names=1)
#rownames(data1)<-colon2hyphen(rownames(data1))
load("./MHL4.RData")
colnames(data)<-gsub(".sorted.clipped.bam.hapInfo.txt.hap","",colnames(data))
colnames(data)<-gsub(".hapInfo.txt.hap","",colnames(data))
saminfo<-read.table("/media/NAS1/shg047/monod/hapinfo/N37Salk.saminfo",sep="\t")

load("gsirlt.RData")              # selected biomarker for each reference
nrow(gsirlt)
gsirlt1<-subset(gsirlt,refMax>0.3 & GSI>0.6)
gsirlt1$group<-as.character(gsirlt1$group)
nrow(gsirlt1)
bio<-gsi2bio(gsirlt1)

Data=data[,grep("STL|N37|age|ENC|new|centenarian|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
Data<-rename(Data)
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
colnames(Data)<-unlist(lapply(colnames(Data),function(x) unlist(strsplit(x,"[.|-]"))[1]))


input<-data.frame(t(apply(Data,1,function(x) tapply(x,colnames(Data),function(x) mean(x,na.rm=T)))))

# for colon cancer plasma
par(mfrow=c(4,1),mar=c(2,5,3,5))
zz<-subset(input, CCP >0.01 & WBC<0.01)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="CCP >0.01 & WBC<0.01")

zz<-subset(input, CCP >0.01 & WBC<0.01 & CCT>0.3)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="CCP >0.01 & WBC<0.01 & CCT>0.3")

zz<-subset(input, CCP >0.01 & WBC<0.01 & Colon>0.3)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="CCP >0.01 & WBC<0.01 & Colon>0.3")

zz<-subset(input, CCP >0.01 & WBC<0.01 & Colon>0.3 & CCT>0.3)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="CCP >0.01 & WBC<0.01 & Colon>0.3 & CCT>0.3")

# for lung cancer plasma
par(mfrow=c(4,1),mar=c(2,5,3,5))
zz<-subset(input, LCP >0.01 & WBC<0.01)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="LCP >0.01 & WBC<0.01")

zz<-subset(input, LCP >0.01 & WBC<0.01 & LCT>0.3)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="LCP >0.01 & WBC<0.01 & LCT>0.3")

zz<-subset(input, LCP >0.01 & WBC<0.01 & Lung>0.3)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="LCP >0.01 & WBC<0.01 & Lung>0.3")

zz<-subset(input, LCP >0.01 & WBC<0.01 & Lung>0.3 & LCT>0.3)
nrow(zz)
colon<-bio[na.omit(match(rownames(zz),rownames(bio))),]
xtt<-table(unlist(lapply(colon$group,function(x) unlist(strsplit(x,"[,]")))))
xtt
barplot(xtt,col=1:10,main="LCP >0.01 & WBC<0.01 & Lung>0.3 & LCT>0.3")




setwd("/home/sguo/Downloads")
zz[is.na(zz)]<-0
par(las=2)
barplot(as.numeric(zz[3,]),names=colnames(zz))
boxplot(zz)

DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
table(colnames(DATA))
gsirlt<-gsi(DATA)
save(gsirlt,file="gsirlt.RData")


nrow(gsirlt)
#gsirlt1<-subset(gsirlt,refMax>0.4 & GSI>0.75)
gsirlt1$group<-as.character(gsirlt1$group)
nrow(gsirlt1)
bio<-gsi2bio(gsirlt1)
accuracy<-ZscorebaseTest(data,bio)
accuracy
#bed<-cor2bed(subset(gsirlt1,group=="Colon")[,1])
#write.table(bed,file="biomarker.bed",col.names=F,row.names=F,sep="\t",quote=F)
#save.image(file="save.RData")

png("newdta-olddata.png")
plot(newdata[,1]~data1[,1],pch=16,cex=0.5,xlab="Dinh's MHL",ylab="Shicheng's MHL")
dev.off()
xx<-data.frame(shicheng=newdata[,1],dinh=data1[,1])
rownames(xx)<-rownames(newdata)
subset(xx,shicheng>0.5 & dinh<0.5)
head(subset(xx,shicheng<0.5 & dinh>0.5))




input<-ccp[match(rownames(bio),rownames(ccp)),]
output<-apply(input,2,function(x) rownames(input)[which(x>0.01)])
rlt1<-lapply(output,function(x) table(unlist(lapply(bio[match(x,rownames(bio)),4],function(x) unlist(strsplit(x,"[,]"))))))
rlt1
input<-ccp[match(rownames(bio),rownames(ccp)),]
output<-apply(input,2,function(x) rownames(input)[which(x>0.01)])
Bio<-bio[match(rownames(input),rownames(bio)),]
back<-table(unlist(lapply(Bio[,4],function(x) unlist(strsplit(x,"[,]")))))
rlt1<-lapply(output,function(x) table(unlist(lapply(bio[match(x,rownames(bio)),4],function(x) unlist(strsplit(x,"[,]"))))))
lapply(rlt1,function(x) x/back)
input<-lcp[match(rownames(bio),rownames(lcp)),]
output<-apply(input,2,function(x) rownames(input)[which(x>0.01)])
rlt2<-lapply(output,function(x) table(unlist(lapply(bio[match(x,rownames(bio)),4],function(x) unlist(strsplit(x,"[,]"))))))
input<-np[match(rownames(bio),rownames(np)),]
output<-apply(input,2,function(x) rownames(input)[which(x>0.01)])
rlt3<-lapply(output,function(x) table(unlist(lapply(bio[match(x,rownames(bio)),4],function(x) unlist(strsplit(x,"[,]"))))))

