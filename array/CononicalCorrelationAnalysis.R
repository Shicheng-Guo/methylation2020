#!/usr/bin/R
setwd("");
install.packages("CCA")
library("CCA")
y<-matrix(rnorm(46858*5791),5790,46857)  # methylation array
x<-matrix(rnorm(5790*200,1,10),5790,200)
dim(x)
dim(y)
res<-rcc(x,y,1,1)

data(nutrimouse)
x=as.matrix(nutrimouse$gene)
y=as.matrix(nutrimouse$lipid)

estim.regul(x,y)
correl=matcor(x,y)
img.matcor(correl,type=1)
img.matcor(correl,type=2)




Chr21_450kMerge.txt.trans

setwd("/home/sguo/methylation")
data<-read.table("Chr21_450kMerge.txt.trans",head=T,row.names=1,sep="\t")
data<-t(data)







