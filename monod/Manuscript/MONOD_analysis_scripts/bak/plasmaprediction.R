#!/usr/bin/env Rscript
library("monod")
library("ggplot2")
library("reshape2")
library("optparse")
# please copy monod package from: /media/Home_Raid1/shg047/monod_1.1.tar.gz and install.packages("monod_1.1.tar.gz")
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input MHL matrix", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

data<-read.table(opt$file,head=T,sep="\t",as.is=T,sep="\t")
# rename the colnames
colnames(data)<-gsub(".sorted.clipped.bam.hapInfo.txt.hap","",colnames(data))
colnames(data)<-gsub(".hapInfo.txt.hap","",colnames(data))
saminfo<-read.table("/media/NAS1/shg047/monod/hapinfo/N37Salk.saminfo",sep="\t")

# Since our plasma RRBS data were low-coverage, only MHL>0.01 regions are reliable to be applied in the analysis (cancer fragment are not samplinga or sequencing randomly)
data<-data[which(unlist(apply(data[,grep(".6P|X6.P|.7P|X7.P|NC.P",colnames(data))],1,function(x) mean(x,na.rm=T)>0.05))),]

# rename the colnames
data<-rename(data)

# GSI for tissue-specific markers
load("gsirlt.RData")              # selected biomarker for each reference

# biomarker further selection procedure
gsirlt1<-subset(gsirlt,refMax>0.3 & GSI>0.6)
gsirlt1$group<-as.character(gsirlt1$group)
bio<-gsi2bio(gsirlt1)
# remove non-WBC marker have higher MHL in WBC
exc<-c()
for(i in 1:nrow(bio)){
  if(length(grep("WBC",bio[i,4]))<1 & ! is.na(as.numeric(unlist(strsplit(as.character(bio[i,7]),","))[10])) & as.numeric(unlist(strsplit(as.character(bio[i,7]),","))[10])>0.01){
    exc<-c(exc,i)
  }
}
bio<-bio[-exc,]

# remove non-WBC marker have higher MHL in WBC
np<-data[,grep("NC.P|NCP|NP",colnames(data))]  
train<-sample(1:75,45)
ncp<-np[,train]
nncp<-np[,(1:75)[-train]]

ccp<-data[,unique(grep(".6P|X6.P|CCP",colnames(data)))]  
lcp<-data[,unique(grep(".7P|X7.P|LCP",colnames(data)))] 

rlt1<-ZscorePredictionPlasma(ccp,ncp,bio,tt=0.01)   # tt is threshold to binary the MHL matrix, ccp: colon cancer plasma, ncp: normal plasma trainning dataset
rlt2<-ZscorePredictionPlasma(lcp,ncp,bio,tt=0.01)   # tt is threshold to binary the MHL matrix, lcp: colon cancer plasma, ncp: normal plasma trainning dataset
rlt3<-ZscorePredictionPlasma(nncp,ncp,bio,tt=0.01)  # tt is threshold to binary the MHL matrix, nncp: normal plasma test dataset, ncp: normal plasma trainning dataset

# Take maximum Z-score as the tissue-of-origin prediction
rlt1$test
rlt2$test
rlt3$test

# how about take P<0.05 as positive detection? 
acc1<-unlist(lapply(apply(rlt1$ccpvalue,2,function(x) rownames(rlt1$ccpvalue)[which(x<0.05)]),function(x) grep("Colon|Intestine",x)[1]))
acc2<-unlist(lapply(apply(rlt2$ccpvalue,2,function(x) rownames(rlt1$ccpvalue)[which(x<0.05)]),function(x) grep("Lung",x)[1]))
acc3<-unlist(lapply(apply(rlt3$ccpvalue,2,function(x) rownames(rlt1$ccpvalue)[which(x<0.05)]),function(x) grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach",x)[1]))
length(na.omit(acc1))/length(acc1)
length(na.omit(acc2))/length(acc2)
length(na.omit(acc3))/length(acc3)
save.image(file="monod.prediction.RData")
