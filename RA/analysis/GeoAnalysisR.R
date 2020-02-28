###############################################
### This file was created to detail the importation and cleaning of the data sets of Rheumatology dataset (SSc, SLE, RA, AS)
### (with associated metadata), and saving them as R-data objects

#### Change this and make it your own!.. preferably should be a folder for data only
setwd("C:/Users/shicheng/Dropbox/Project/methylation/RA")
### You might have to install these packages using two lines:
# source('http://bioconductor.org/biocLite.R')
# biocLite('packageName')
# If they are installed, load them:
library("GEOquery")
library("affy")
library("simpleaffy")
### Get each GSE series from GEO site (this will download the "series matrix file" 
#can take awhile so once downloaded save as .R object for later)

##################################################################
### Get data from GEObase, extract expressions, (do once)
##################################################################
GSE27895 <- getGEO("GSE27895")
show(GSE27895) ## 8 APL and 2 healthy marrow
save(GSE27895, file="GSE27895_matrix.Rdata")
getwd()
### Extract expression matrices (turn into data frames at once) 
load("GSE27895_matrix.Rdata")
dat <- as.data.frame(exprs(GSE27895[[1]]))
### Obtain the meta-data for the samples and rename them perhaps?
meta <- pData(phenoData(GSE27895[[1]]))
meta$characteristics_ch1.1
pheno<-as.numeric(meta[match(colnames(dat),rownames(meta)),]$characteristics_ch1.1)-1
meta$group<-c(rep('Case', 11), rep('Con', 12))

pvalue<-apply(dat,1,function(x) t.test(x[1:11],x[12:23])$p.value)
save(pvalue,file="GSE27895_Pvalue.Rdata")
sort(pvalue)[1:100]
0.05/length(pvalue)
