setwd("C:\\Users\\User\\Dropbox\\Project\\methylation\\esophageal_carcinoma\\phase1\\ZNF132")
data<-read.table("meth-exp.esca.txt",head=T)
summary(lm(data$MF~log(data$EXP)))
plot(log(data$EXP) ~ MF, data = data,col=as.numeric(as.factor(unlist(lapply(as.character(data$Target),function(x) substr(x,1,1)))))+1,xlab="Average Methylation Level")
abline(lm(log(data$EXP)~data$MF), col="blue",lwd=2) # regression line (y~x) 
legend("topright",legend=c("Normal","Cancer"),col=c(2,3),pch=1,bty="n")
           