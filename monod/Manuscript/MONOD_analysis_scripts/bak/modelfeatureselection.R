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
Data=data[,grep("STL|N37|age|ENC|new|centenarian|CTT|HCT|X7.T|X6.T",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
Data<-rename(Data)
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
colnames(Data)<-unlist(lapply(colnames(Data),function(x) unlist(strsplit(x,"[.|-]"))[1]))
gsirlt<-gsi(Data)
save(gsirlt,"gsirlt.RData")        # save biomarker for further analysis
