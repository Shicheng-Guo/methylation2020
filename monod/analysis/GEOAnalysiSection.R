#####################################################################
###   TiTle:  Meta-Analysis of GEO Dataset to reveal now acknowledge
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Section 1. Methylation Dataset Collection 
###   Section 2. PCA Analysis
###   Section 3. Differential Analysis
###   Section 4. Pathway Analysis
#####################################################################

####################################################################################################################################
### Section 1. Methylation Dataset Collection: Get data from GEObase, extract expressions, (do once)
####################################################################################################################################
setwd('/home/sguo/monod/data/geo')
### You might have to install these packages using two lines:
##source('http://bioconductor.org/biocLite.R')
##biocLite('simpleaffy')
# If they are installed, load them:
library("GEOquery")
library("affy")
library("simpleaffy")
### Get each GSE series from GEO site (this will download the "series matrix file" 
#can take awhile so once downloaded save as .R object for later)

######################
### Get data from GEObase, extract expressions, (do once)
######################
GSE56044 <- getGEO("GSE56044",destdir="/home/sguo/monod/data/geo")
save(GSE56044, file="GSE56044_matrix.Rdata")
GSE39279 <- getGEO("GSE39279",destdir="/home/sguo/monod/data/geo")
save(GSE39279, file="GSE39279_matrix.Rdata")
GSE52401 <- getGEO("GSE52401",destdir="/home/sguo/monod/data/geo")
save(GSE52401, file="GSE52401_matrix.Rdata")
GSE42752 <- getGEO("GSE42752",destdir="/home/sguo/monod/data/geo")
save(GSE42752, file="GSE42752_matrix.Rdata")
GSE27895 <- getGEO("GSE27895",destdir="/home/sguo/monod/data/geo")
save(GSE27895, file="GSE27895_matrix.Rdata")
GSE37066 <- getGEO("GSE37066",destdir="/home/sguo/monod/data/geo")
save(GSE37066, file="GSE37066_matrix.Rdata")
GSE61461 <- getGEO("GSE61461",destdir="/home/sguo/monod/data/geo")
save(GSE61461, file="GSE61461_matrix.Rdata")
GSE50192 <- getGEO("GSE50192",destdir="/home/sguo/monod/data/geo")
save(GSE50192, file="GSE50192_matrix.Rdata")
GSE48472 <- getGEO("GSE48472",destdir="/home/sguo/monod/data/geo")
save(GSE48472, file="GSE48472_matrix.Rdata")
GSE64491 <- getGEO("GSE64491",destdir="/home/sguo/monod/data/geo")
save(GSE64491, file="GSE64491_matrix.Rdata")


######################
### load the RData
######################

load("GSE56044_matrix.Rdata")
load("GSE39279_matrix.Rdata")
load("GSE52401_matrix.Rdata")
load("GSE42752_matrix.Rdata")
load("GSE27895_matrix.Rdata")

### Extract expression matrices (turn into data frames at once) 
dat_GSE56044 <- as.data.frame(exprs(GSE56044[[1]]))
dat_GSE39279 <- as.data.frame(exprs(GSE39279[[1]]))
dat_GSE52401 <- as.data.frame(exprs(GSE52401[[1]]))
dat_GSE42752 <- as.data.frame(exprs(GSE42752[[1]]))
dat_GSE27895 <- as.data.frame(exprs(GSE27895[[1]]))

### Obtain the meta-data for the samples and rename them perhaps?
p_GSE56044 <- pData(phenoData(GSE56044[[1]]))
p_GSE39279 <- pData(phenoData(GSE39279[[1]]))
p_GSE52401 <- pData(phenoData(GSE52401[[1]]))
p_GSE42752 <- pData(phenoData(GSE42752[[1]]))
p_GSE27895 <- pData(phenoData(GSE27895[[1]]))

p_GSE56044$source_name_ch1
p_GSE39279$source_name_ch1
p_GSE52401$characteristics_ch1.2
p_GSE42752$source_name_ch1
p_GSE27895$characteristics_ch1.1


table(p_GSE56044$source_name_ch1)       # GSE56044: 124 lung cancer and 12 normal tissues
table(p_GSE39279$source_name_ch1)       # GSE39279: 322 adenocarcinoma and 122 squamous lung cancer
table(p_GSE52401$characteristics_ch1.2) # GSE52401: 244 adjacent lung tissues
table(p_GSE42752$source_name_ch1)       # GSE42752:  22 colon  adenocarcinoma, 22 corresponding normal and 19 cancer-unrelated normal colon tissue
table(p_GSE27895$characteristics_ch1.1) # GSE42752:  11 RA and 12 Normal CD4+

#####################
### Rename samples
#####################
dat_GSE56044
dat_GSE39279
dat_GSE52401
dat_GSE42752
dat_GSE27895

#########################
### phenotype regroup
#########################
PhenRecode<-function(dat,p,phen="source_name_ch1"){
nphen<-match(phen,colnames(p))
phen<-as.numeric(p[match(colnames(dat),p$geo_accession),nphen])
phen
}

phen_GSE56044<-PhenRecode(dat_GSE56044,p_GSE56044,phen="source_name_ch1")
phen_GSE39279<-PhenRecode(dat_GSE39279,p_GSE39279,phen="source_name_ch1")
phen_GSE52401<-PhenRecode(dat_GSE52401,p_GSE52401,phen="characteristics_ch1.2")
phen_GSE42752<-PhenRecode(dat_GSE42752,p_GSE42752,phen="source_name_ch1")
phen_GSE27895<-PhenRecode(dat_GSE27895,p_GSE27895,phen="characteristics_ch1.1")

# group 1: cancer
lungcancer_GSE56044<-dat_GSE56044[,phen_GSE56044="1"] 
lungcancer_GSE39279<dat_GSE39279[,phen_GSE39279="1"]  


# group 2: adjacent normal
lungadjnormal_GSE56044<-dat_GSE56044[,phen_GSE56044=="2"] 
lungadjnormal_GSE52401<-dat_GSE52401[,phen_GSE52401=="1"]


# group 3: abslote normal
colonabsnormal_GSE42752<-dat_GSE42752[,phen_GSE42752="2"]
cd4normal_GSE27895<-dat_GSE27895[,phen_GSE27895=="1"]

# group 4: RA
cd4rheumatoid_GSE27895<-dat_GSE27895[,phen_GSE27895=="2"] 

# group 5: SLE


# group 5: Gout


####################################################################################################################################
### Section 2. Differential Methylation Analysis
####################################################################################################################################

# Define the function to test the DMR
DMRtest<-function(case,control){
data<-data.matrix(data.frame(case,control))
x1<-ncol(case)
x2<-ncol(control)
x3<-ncol(data)
# wilcox test and median, beta, fold fold 
pvalue<-apply(data,1,function(x) p<-wilcox.test(x[1:x1],x[(x1+1):x3],na.omit=T)$p.value)
med1<-apply(data,1,function(x) median(x[1:580],na.rm=T))
med2<-apply(data,1,function(x) median(x[581:960],na.rm=T))
beta<-med1-med2
ratio<-med1/med2
hypercase<-apply(data,1,function(x) sum(x[1:580]>0.3,na.rm=T)/length(na.omit(x)))
hypercon<-apply(data,1,function(x) sum(x[581:960]>0.3,na.rm=T)/length(na.omit(x)))
output<-data.frame(pvalue,med1,med2,beta,ratio,hypercase,hypercon)
sig<-subset(output,abs(beta)>0.1 & ratio>1.2 & hypercase>0.6)
write.table(output,file="lungcance.methylation.rlt.txt",col.names=NA,row.names=T,quote=F,sep="\t")
write.table(sig,file="lungcance.sig.dmr.rlt.txt",col.names=NA,row.names=T,quote=F,sep="\t")

return(output)
}



# Lung Cancer Project
lungcancer<-data.frame(lungcancer_GSE56044,lungcancer_GSE39279)
lungnormal<-data.frame(lungadjnormal_GSE56044,lungadjnormal_GSE52401)

# Rheumatoid project
cd4RA<-data.frame(cd4rheumatoid_GSE27895)
cd4normal<-data.frame(cd4normal_GSE27895)

# Pancreatic Cancer Project


# 



####################################################################################################################################
### Section 3. Principle Componment Analysis
####################################################################################################################################

PCAPlot<-function(data,pheno,output,multifigure=T){
  pca <- prcomp(data,center=T,scale = F)  # Here, input file: row is individual and column is variable
  outputfile=paste(output,".pdf",sep="")
  pdf(outputfile)
  if(multifigure){
    par(mfrow=c(2,2),mar=c(4,4,4,4))  
  }
  plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
  axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
  var<-c()
  for(i in 1:length(pca$sdev)){var[i]<-sum((pca$sdev[1:i])^2)/sum((pca$sdev)^2)}
  plot(var,ylab="total variance",xlab="number of principle components",lwd=2,type="l")
  abline(h=0.8,col="grey",lty=2)
  abline(v=which(var>0.8)[1],col="grey",lty=2)
  scores <- data.frame(pheno, pca$x[,1:3])
  col = as.numeric(as.factor(pheno))
  plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.8,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n")
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n")
  dev.off()
}




