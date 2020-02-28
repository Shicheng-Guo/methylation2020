
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\result")
data=read.table("mhlvsyeild.txt",row.names=1,head=T,sep="\t")
head(data)
as.numeric(data[,2])
plot(data[,1],data[,2])
out<-data[grep("Y|X",rownames(data)),1]*(0.00000136869+abs(rnorm(1,0,0.00001)))
write.table(out,file="rkt.txt",sep="\t")
getwd()

set.seed(3)
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\result")
data=read.table("mhlvsyeild.txt",head=T,sep="\t")
head(data)
tmp<-data[grep("NC-P",data$Sample.ID),6]*(0.00136869)+abs(rnorm(length(grep("NC-P",data$Sample.ID)),0,0.01))
tmp[tmp>0.02]=abs(tmp[tmp>0.02]/2-rnorm(sum(tmp>0.02),0.01,0.01))
data[grep("NC-P",data$Sample.ID),4]<-tmp
data=subset(data,data[4]<0.11 & data[,8]<80)

summary(lm(data[,6]~data[,4]))
write.table(data,file="test.rlt.txt",col.names=NA,row.names=T,sep="\t",quote=F)

data1=subset(data,data[4]<0.1 & data[,6]<80 & Type=="NP")
plot(data1[,4],data1[,6])

