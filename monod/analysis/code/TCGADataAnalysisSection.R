

setwd("/home/sguo/monod/data")
load("TCGA.COAD.450K.RData")
load("TCGA.LUAD.450K.RData")
load("TCGA.LUSC.450K.RData")
load("TCGA.PAAD.450K.RData")

dim(coad)
dim(luad)
dim(lusc)
dim(paad)

table(substr(colnames(coad),14,15))
table(substr(colnames(luad),14,15))
table(substr(colnames(lusc),14,15))
table(substr(colnames(paad),14,15))

phen1<-(substr(colnames(coad),14,15))
phen2<-(substr(colnames(luad),14,15))
phen3<-(substr(colnames(lusc),14,15))
phen4<-(substr(colnames(paad),14,15))

type1<-c(which(phen1=="01"), which(phen1=="11"))
type2<-c(which(phen2=="01"), which(phen2=="11"))
type3<-c(which(phen3=="01"), which(phen3=="11"))
type4<-c(which(phen4=="01"), which(phen4=="11"))

coad<-coad[,type1]
luad<-luad[,type2]
lusc<-lusc[,type3]
paad<-paad[,type4]


