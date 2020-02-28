# for ips methylatin 450K analysis
library("GEOquery")
GSE76269 <- getGEO("GSE76269")
data <- as.data.frame(exprs(GSE76269[[1]]))
phen <- pData(phenoData(GSE76269[[1]]))

data[1:5,1:5]
phen[1:5,1:5]

normal=phen[phen$`health state:ch1`=="normal",]
cancer=phen[phen$`health state:ch1`=="cancer",]

phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.7),"[:]"),function(x) as.numeric(unlist(x)[2]))  # status 1:control, 2:scz
phen1[phen1==1]<-"Normal"
phen1[phen1==2]<-"schizophrenia"


phen2<-sapply(strsplit(as.character(phen$characteristics_ch1),"[:]"),function(x) (unlist(x)[2]))  # gender

data1=na.omit(data)

PCAPlot(t(data1),phen1,output="GSE41169.scz.normal.pdf",multifigure=T)  # status
PCAPlot(t(data1),phen2,output="GSE41169.gender.pdf",multifigure=T)  # gender

newphen=data.frame(phen$title,phen$`health state:ch1`)
newphen
