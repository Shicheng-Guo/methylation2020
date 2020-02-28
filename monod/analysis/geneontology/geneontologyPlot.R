
pdf("Gene.ontology.pdf")

par(mfrow=c(4,2))
data<-read.delim("C:\\Users\\shicheng\\Downloads\\save.txt",header=FALSE)
data
data[,2]<--log(data[,2],10)
data
data<-data[order(data[,2]),]
data
par(las=1)
par(mar=c(4,20,2,2)) # increase y-axis margin.
barplot(data[,2],horiz=TRUE,names.arg=data[,1],cex.names=0.8,col="blue",xlim=c(0,3),xlab="z")


data<-read.delim("C:\\Users\\shicheng\\Downloads\\save2.txt",header=FALSE)
data
data[,2]<--log(data[,2],10)
data
data<-data[order(data[,2]),]
data
par(las=1)
par(mar=c(2,20,2,2)) # increase y-axis margin.
barplot(data[,2],horiz=TRUE,names.arg=data[,1],cex.names=0.8,col="blue",xlim=c(0,3),xlab="-log10(Benjamini P value)")

dev.off()


file=list.files(path="/home/shg047/monod/hap/wgbs/All_chromosomes_combined/",pattern="*hapInfo.txt.pairwiseR2")

data<-c()
for(i in 1:length(file)){

  tmp<-read.table(file[i],skip=1)
  data<-data.frame(data,tmp[,3])
  
  
}



