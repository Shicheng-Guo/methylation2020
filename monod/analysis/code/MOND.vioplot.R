setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\MHL-March30-2016")
data<-read.table("Figure4-lungcancer.data.txt",head=T,row.names=1,sep="\t",check.names=F)

# install.packages("vioplot")
library("vioplot")
par(mfrow=c(2,1))
normal.plasma<-data[,grep("NP-Kun",colnames(data))]
cancer.plasma<-data[,grep("LCP",colnames(data))]
colnames(data)
vioplot(normal.plasma,cancer.plasma)



mu<-2
si<-0.6
bimodal<-c(rnorm(1000,-mu,si),rnorm(1000,mu,si)) 
uniform<-runif(2000,-4,4)
normal<-rnorm(2000,0,3)
vioplot(bimodal,uniform,normal)
boxplot(bimodal,uniform,normal)

# add to an existing plot
x <- rnorm(100)
y <- rnorm(100)
plot(x, y, xlim=c(-5,5), ylim=c(-5,5))
vioplot(x, col="tomato", horizontal=TRUE, at=-4, add=TRUE,lty=2, rectCol="gray")
vioplot(y, col="cyan", horizontal=FALSE, at=-4, add=TRUE,lty=2)
