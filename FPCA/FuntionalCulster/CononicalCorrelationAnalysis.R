#!/usr/bin/R
setwd("");
install.packages("CCA")
library("CCA")
y<-matrix(rnorm(20000),200,100)
x<-matrix(rnorm(40000,1,10),200,200)
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


