setwd("/home/sguo/dyh/idat")
map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",sep="\t",as.is=T)
sam<-read.table("saminfo.txt",head=T,sep="\t",as.is=T)
sam2<-read.table("saminfo2.txt",head=T,sep="\t",as.is=T)

dms<-read.table("dms.txt",sep="\t",as.is=T,head=T)
load("myNorm.RData")

ssid<-sam[match(colnames(myNorm$beta)[1:12],sam$Sample_Name),]$SSID

x1<-sam2[match(ssid,sam2$SSID)]$History
x2<-sam2[match(ssid,sam2$SSID)]$RF
x3<-sam2[match(ssid,sam2$SSID)]$ESR
x4<-sam2[match(ssid,sam2$SSID)]$GJZZS
x5<-sam2[match(ssid,sam2$SSID)]$GJYTS
x6<-sam2[match(ssid,sam2$SSID)]$HZPG
x7<-sam2[match(ssid,sam2$SSID)]$DAS28

newdata<-myNorm$beta[,1:12]
newdata<-newdata[na.omit(match(dms[,1],rownames(newdata))),]

p1<-apply(newdata,1,function(x) summary(lm(x1~x))$coefficients[2,4])
p2<-apply(newdata,1,function(x) summary(lm(x2~x))$coefficients[2,4])
p3<-apply(newdata,1,function(x) summary(lm(x3~x))$coefficients[2,4])
p4<-apply(newdata,1,function(x) summary(lm(x4~x))$coefficients[2,4])
p5<-apply(newdata,1,function(x) summary(lm(x5~x))$coefficients[2,4])
p6<-apply(newdata,1,function(x) summary(lm(x6~x))$coefficients[2,4])
p7<-apply(newdata,1,function(x) summary(lm(x7~x))$coefficients[2,4])

gene<-map[match(rownames(newdata),map[,4]),5]

rlt1<-cbind(newdata,gene,p1)
rlt2<-cbind(newdata,gene,p2)
rlt3<-cbind(newdata,gene,p3)
rlt4<-cbind(newdata,gene,p4)
rlt5<-cbind(newdata,gene,p5)
rlt6<-cbind(newdata,gene,p6)
rlt7<-cbind(newdata,gene,p7)

subset(rlt1,p1<0.001)
subset(rlt2,p2<0.001)
subset(rlt6,p6<0.001)

subset(rlt3,p3<0.005)
subset(rlt4,p4<0.005)
subset(rlt5,p5<0.005)
subset(rlt7,p7<0.005)



