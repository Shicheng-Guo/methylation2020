# 2017-05-10

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("impute")  
#install.packages("gplots")
#install.packages("RColorBrewer")
#install.packages("grDevices")
library("gplots")
library("RColorBrewer")
library("grDevices")
library("impute")


gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  GSI<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    GSI<-c(GSI,gsit)
    #print(paste(gmax,gsit,sep=" "))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=GSI)
  return(rlt)
}

TopGSIByCategory<-function(gsi,top=150){
  GSIRlt<-c()
  group<-names(table(gsi$group))
  rank<-c(rep(top,length(group)))
  for (i in 1:length(group)){
    subset=gsi[which(gsi$group==group[i]),]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
  }
  return(na.omit(GSIRlt))
}

setwd("/home/shg047/oasis/monod/hapinfo")
saminfo<-read.table("/home/shg047/oasis/monod/saminfo/N37Salk.saminfo",sep="\t")
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
# feature reduction (WGBS+RRBS missing<60, each plasma category(colon,lung and normal)<60%, Plasma missing<50%)
rm1<-which(apply(data[,grep("X7.P|X6.P|.6P|.7P|NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.5))
rm2<-which(apply(data[,grep("X6.P|.6P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
rm3<-which(apply(data[,grep("X7.P|.7P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
rm4<-which(apply(data[,grep("NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
rm5<-which(apply(data,1,function(x) sum(is.na(x))/length(x)>0.6))
rm<-unique(c(rm1,rm2,rm3,rm4,rm5))
data<-data[-rm,]

checkdata<-data[match(rownames(bio),rownames(data)),]

# copy our biomarker form supplementary table to: /home/sguo/Dropbox/Project/methylation/monod/Manuscript/MONOD_analysis_scripts/biomarker2.txt
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)<-gsub("[.]","-",colnames(Data))
colnames(Data)[grep("age|new|centenarian",colnames(Data))]<-"WBC"
colnames(Data)[grep("X7.T|X6.T|SRX|CTT",colnames(Data))]<-"CT"
colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
# for tissue-specific biomarkers
library("impute")
library("preprocessCore")
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WB",colnames(Data))] 
DATA<-Data[match(rownames(bio),rownames(Data)),grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WB",colnames(Data))] 
colnames(DATA)<-colnames(Data)[grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WB",colnames(Data))]
for(i in names(table(colnames(DATA)))){
  DATA[,colnames(DATA)==i]<-normalize.quantiles(data.matrix(DATA[,colnames(DATA)==i]))
}
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
gsirlt<-gsi(DATA)
topgsi<-TopGSIByCategory(gsirlt,top=300)   


# variable selection
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
Bio<-bio[sample(1:nrow(bio),0.9*nrow(bio)),]  # apply parts of the tissue-specific biomarkers so that we can evalate the stability
DATA<-Data[match(rownames(bio),rownames(Data)),grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WB",colnames(Data))] 
matchrlt<-apply(DATA,2,function(x) names(table(bio$V5))[which.max(tapply(x,bio$V5,function(x) sum(x>0.5,na.rm=T)))])
for(i in names(table(bio$V5))){
  acc=sum(matchrlt[grep(i,names(matchrlt))]==i)/length(grep(i,names(matchrlt)))
  print(c(i,acc))
}

## Result: Eventually, You will get the accuray as the following:
# [1] "Brain" "1"    
# [1] "Colon" "1"    
# [1] "Intestine" "0.75"     
# [1] "Kidney" "0.5"   
# [1] "Liver" "1"    
# [1] "Lung" "0.8" 
# [1] "Pancreas" "1"       
# [1] "Spleen"            "0.666666666666667"
# [1] "Stomach" "0.5"    
# [1] "WBC" "NaN"