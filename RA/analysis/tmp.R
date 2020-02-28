

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\RA\\analysis\\paperOARA\\")
d1<-read.table("DMG.txt")
head(d1)

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\RA\\analysis\\")

d2<-read.table("hyper.gene.txt")
d3<-read.table("hypo.gene.txt")

d2<-rbind(d2,d3)


length(na.omit(match(d1[,1],d2[,1])))
