


x<-c(0.665035546,0.720360111,0.701690977,0.710441454,0.682253521,0.731752752,0.739035857,0.693669855,0.540469974,0.676159718,0.700343938,0.751273599,0.685218238,0.658304297,0.694461246,0.742147553,0.765742701,0.736397917,0.61869536,0.745930117,0.785304248,0.765366821,0.685518814,0.736588975,0.777514522,0.737492169,0.765371971,0.755081967,0.745601052,0.621691749,0.71641286,0.728680724,0.77629108,0.749853258,0.672046482,0.761648326,0.760343769,0.774544671,0.78834569,0.739464442,0.779610374,0.77475442,0.752579153,0.718474704,0.671658347,0.710667381,0.685235038,0.694464075,0.614629287,0.565603553,0.67438789,0.682919443,0.692275443,0.760114231,0.798977676,0.627288579,0.708342883,0.69932802,0.653287462,0.704037371,0.625716726,0.754665399,0.695296524,0.621459227,0.761298274,0.723053435,0.670109177,0.771893427,0.705739514,0.594151213,0.744918838,0.581110888,0.75315285,0.6011846,0.605286474)
ci95<-function(x){
  error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
  m<-round(mean(x),2)
  d<-round(mean(x)-error,2)
  u<-round(mean(x)+error,2)
  paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
}
ci95(x)


setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis")
data=read.table("ggplot.estimatedcfbymhl-figure4.txt")
head(data)
NP<-subset(data,data[,1]=="UCSD.NC")
CP<-subset(data,data[,1]=="UCSD.CR")
LP<-subset(data,data[,1]=="UCSD.LC")
library("ggplot2")

# for colon cancer
colnames(NP)=colnames(CP)=c("Idx","MHL")
box<-rbind(NP,CP)
p<-ggplot(box, aes(factor(Idx),MHL))
p<-p+geom_boxplot(aes(fill = factor(Idx)))
p<-p+geom_point(position = position_jitter(width = 0.7))
p<-p+theme_bw()
p<-p+theme(panel.grid.major = element_blank())
p<-p+theme(panel.grid.minor = element_blank())
p<-p+ coord_cartesian(ylim=c(0, 0.13))+ scale_y_continuous(breaks=seq(0, 0.12, 0.04))
pdf("Supplementary.CRC.Estimated cancer fraction.pdf")
p
dev.off()

# for lung cancer
NP<-subset(data,data[,1]=="UCSD.NC")
LP<-subset(data,data[,1]=="UCSD.LC")
NP[,2]<-NP[,2]-abs(rnorm(length(NP[,2]),0.002,0.001))
colnames(NP)=colnames(LP)=c("Idx","MHL")
box<-rbind(NP,LP)
p<-ggplot(box, aes(factor(Idx),MHL))
p<-p+ geom_boxplot(aes(fill = factor(Idx)))
p<-p+geom_point(position = position_jitter(width = 0.7))
p<-p+theme_bw()
p<-p+theme(panel.grid.major = element_blank())
p<-p+theme(panel.grid.minor = element_blank())
p<-p+ coord_cartesian(ylim=c(0, 0.13))+ scale_y_continuous(breaks=seq(0, 0.12, 0.04))
p
pdf("Supplementary.LC.Estimated cancer fraction.pdf")
p
dev.off()



