

### Phenotype

setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\esophageal_carcinoma\\ZNF132")
phen<-read.table("ESCA.phenotype.txt",head=T)
head(phen)
T<-unlist(lapply(phen$TNM,function(x) substr(x,1,2)))
N<-unlist(lapply(phen$TNM,function(x) substr(x,3,4)))
M<-unlist(lapply(phen$TNM,function(x) substr(x,5,6)))
output<-data.frame(phen,T,N,M)
head(output)
write.table(output,file="ESCA.phenotype.new.txt",quote=F,sep="\t")

### DNA methylation
setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\esophageal_carcinoma\\ZNF132")
data<-read.table("meth-exp.esca.txt",head=T,sep="\t",row.names = 1)
head(data)
summary(lm(data$MF~log(data$EXP)))
plot(log(data$EXP) ~ MF, data = data,col=as.numeric(as.factor(unlist(lapply(as.character(data$Target),function(x) substr(x,1,1)))))+1,xlab="Average Methylation Level")
abline(lm(log(data$EXP)~data$MF), col="blue",lwd=2) # regression line (y~x) 
legend("topright",legend=c("Normal","Cancer"),col=c(2,3),pch=1,bty="n")
