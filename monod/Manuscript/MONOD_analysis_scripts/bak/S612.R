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

CvSampling<- function(Nobs=29,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

topgsi2bio<-function(topgsi){
  cor2bed<-function(cor){
    cor<-as.character(cor)
    a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
    bed<-matrix(a,ncol=3,byrow=T)
    return(data.frame(bed))
  }
  bio<-data.frame(cor2bed(topgsi[,1]),topgsi[,2:3])
  rownames(bio)<-topgsi[,1]
  colnames(bio)<-paste("V",2:6,sep="")
  return(bio)
}

rename<-function(data){
  Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P",colnames(data))]
  colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
  colnames(Data)<-gsub("[.]","-",colnames(Data))
  colnames(Data)[grep("age|new|centenarian",colnames(Data))]<-"WBC"
  colnames(Data)[grep("X7.T|X6.T|SRX|CTT",colnames(Data))]<-"CT"
  colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
  colnames(Data)[grep("X6.P|RRBS.6P",colnames(Data))]<-"CCP"
  colnames(Data)[grep("X7.P|RRBS.7P",colnames(Data))]<-"LCP"
  return(Data)
}

frm<-function(data){
  # feature reduction (WGBS,RRBS and each caterogy missing<60,Plasma missing<50%,low variation removed)
  # features which are missing in whole reference samples will be omit since this will caused miss-classification (N<=4,data qualitye dependent)
  rm1<-which(apply(data[,grep("X7.P|X6.P|.6P|.7P|NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.5))
  rm2<-which(apply(data[,grep("X6.P|.6P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
  rm3<-which(apply(data[,grep("X7.P|.7P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
  rm4<-which(apply(data[,grep("NC.P",colnames(data))],1,function(x) sum(is.na(x))/length(x)>0.6))
  rm5<-which(apply(data,1,function(x) sum(is.na(x))/length(x)>0.6))
  rm<-unique(c(rm1,rm2,rm3,rm4,rm5))
  return(rm)
}

setwd("/home/shg047/oasis/monod/hapinfo")
data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
saminfo<-read.table("/home/shg047/oasis/monod/saminfo/N37Salk.saminfo",sep="\t")
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
#load("/oasis/tscc/scratch/shg047/monod/hapinfo/Depth4.RData")
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 

rm<-frm(data)
data<-data[-rm,]
#depth<-depth[-rm,]
# copy our biomarker form supplementary table to: /home/sguo/Dropbox/Project/methylation/monod/Manuscript/MONOD_analysis_scripts/biomarker2.txt (download from supplementary)
Data<-rename(data)
colon<-subset(bio,V5=="Colon")
DATA<-data[match(rownames(colon),rownames(data)),]
DATA<-rename(DATA)
DATA<-DATA[,colnames(DATA)=="CCP"]
apply(DATA,1,function(x) sum(x>0.1,na.rm=T))
TDTD<-data[match(names(sort(apply(DATA,1,function(x) sum(x>0.3,na.rm=T)),decreasing=T)[1:10]),rownames(data)),]
TDTD<-rename(TDTD)

# for tissue-specific biomarkers
#library("impute")
#library("preprocessCore")
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
DATA<-DATA[match(rownames(bio),rownames(DATA)),]
colnames(DATA)<-colnames(Data)[grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))]
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))

gsirlt<-gsi(DATA)
for(i in names(table(colnames(DATA)))){
  
}
apply(DATA,1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T)))
topgsi<-TopGSIByCategory(gsirlt,top=500)   

bio<-topgsi2bio(topgsi)
Bio<-bio[sample(1:nrow(bio),0.9*nrow(bio)),]  # apply parts of the tissue-specific biomarkers so that we can evalate the stability
DATA<-Data[match(rownames(bio),rownames(Data)),grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
matchrlt<-apply(DATA,2,function(x) names(table(bio$V5))[which.max(tapply(x,bio$V5,function(x) sum(x>0.1,na.rm=T)))])
for(i in names(table(bio$V5))){
  acc=sum(matchrlt[grep(i,names(matchrlt))]==i)/length(grep(i,names(matchrlt)))
  print(c(i,acc))
}
#[1] "Brain" "1"
#[1] "Colon" "1"
#[1] "Intestine" "1"
#[1] "Kidney" "1"
#[1] "Liver" "1"
#[1] "Lung" "0.8"
#[1] "Pancreas" "1"
#[1] "Spleen" "1"
#[1] "Stomach""0.666666666666667"
#[1] "WBC" "1"

#### Figure 5C
set.seed(2)
out1<-grep(".6P|X6.P",colnames(data))[sample(1:30,5)]
out2<-grep(".7P|X7.P",colnames(data))[sample(1:30,5)]
out3<-grep("NC.P",colnames(data))[sample(1:75,5)]
test<-data[,c(out1,out2,out3)]
test<-test[match(rownames(bio),rownames(test)),]
# Cross-validation process
Data<-data[,-c(out1,out2,out3)]
# automatically select best threshold with 5-fold cross-validation for colon plasma/lung cancer/normal plasma together.
input<-Data[match(rownames(bio),rownames(Data)),]
set.seed(51)
k=2                    # split the sample to two parts, one is for train and the other is for test. 
acc1<-c()
acc2<-c()
acc3<-c()
Best<-c()
Samping<-CvSampling(24,k)
Lnum1<-c()
Lnum2<-c()
Lnum3<-c()
for(i in 1:k){
  Num<-c()
  for(j in seq(0,1,0.01)){
    counts1<-apply(input[,grep(".6P|X6.P",colnames(Data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    counts2<-apply(input[,grep(".7P|X7.P",colnames(Data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    counts3<-apply(input[,grep("NC.P",colnames(Data))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
    num<-data.frame(id=j,
                    c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                    c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                    c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  best<-Num[which.max(rowSums(Num[,2:4])),1]
  countm1<-apply(input[,grep(".6P|X6.P",colnames(Data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm2<-apply(input[,grep(".7P|X7.P",colnames(Data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  countm3<-apply(input[,grep("NC.P",colnames(Data))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
  Lnum1<-c(Lnum1,c=sum(apply(countm1,2,function(x) which.max(x)==2)))
  Lnum2<-c(Lnum2,c=sum(apply(countm2,2,function(x) which.max(x)==6)))
  Lnum3<-c(Lnum3,c=sum(apply(countm3,2,function(x) which.max(x)==10)))
  
  acc1<-rbind(acc1,c(best,sum(apply(countm1,2,function(x) which.max(x)==2))/(length(Samping[[i]]$test))))
  acc2<-rbind(acc2,c(best,sum(apply(countm2,2,function(x) which.max(x)==6))/(length(Samping[[i]]$test))))
  acc3<-rbind(acc3,c(best,sum(apply(countm3,2,function(x) which.max(x)==10))/(length(Samping[[i]]$test))))
  
  Best<-c(Best,best)
}
Best
cc1<-apply(test[,grep(".6P|X6.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
cc2<-apply(test[,grep(".7P|X7.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
cc3<-apply(test[,grep("NC.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
ccc1=sum(apply(cc1,2,function(x) which.max(x)==2))/5
ccc2=sum(apply(cc2,2,function(x) which.max(x)==6))/5
ccc3=sum(apply(cc3,2,function(x) which.max(x)==10))/5
ccc1
ccc2
ccc3

######################################################################
################  Model Stability ####################################
######################################################################
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
starlt<-stability(data,bio,rep=200)
save(starlt,file="starlt.RData")
# starlt
#        [,1] [,2] [,3] [,4]
# rlt 0.3605  0.8  0.8  0.8
# rlt 0.3390  0.8  1.0  1.0

senthresplot<-function(starlt){
  acc<-data.frame(x=acc[,1],y=acc[,2])
  acc$bin<- cut(acc[,1], c(seq(0,1,0.01)))
  ggplot(acc) + geom_boxplot(aes(bin, y))
}

acc=starlt[,c(1,2)]
senthresplot(acc)
ggsave("colon-threshold-acc.pdf")

acc=starlt[,c(1,3)]
senthresplot(acc)
ggsave("lung-threshold-acc.pdf")

acc=starlt[,c(1,4)]
senthresplot(acc)
ggsave("normal-threshold-acc.pdf")
stability<-function(data,bio,rep){
  Rlt<-c()
  for(loop in 1:rep){
    out1<-grep(".6P|X6.P",colnames(data))[sample(1:30,5)]
    out2<-grep(".7P|X7.P",colnames(data))[sample(1:29,5)]
    out3<-grep("NC.P",colnames(data))[sample(1:75,5)]
    test<-data[,c(out1,out2,out3)]
    test<-test[match(rownames(bio),rownames(test)),]
    Data<-data[,-c(out1,out2,out3)]
    # automatically select best threshold with 5-fold cross-validation for colon plasma/lung cancer/normal plasma together.
    input<-Data[match(rownames(bio),rownames(Data)),]
    k=2
    acc1<-c()
    acc2<-c()
    acc3<-c()
    Best<-c()
    Samping<-CvSampling(24,k)
    Lnum1<-c()
    Lnum2<-c()
    Lnum3<-c()
    for(i in 1:k){
      Num<-c()
      for(j in seq(0,1,0.001)){
        counts1<-apply(input[,grep(".6P|X6.P",colnames(input))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
        counts2<-apply(input[,grep(".7P|X7.P",colnames(input))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
        counts3<-apply(input[,grep("NC.P",colnames(input))[Samping[[i]]$train]],2,function(x) tapply(x,bio$V5,function(x) sum(x>j,na.rm=T)))
        num<-data.frame(id=j,
                        c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                        c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                        c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
        Num<-rbind(Num,num)
      }
      best<-Num[which.max(rowSums(Num[,2:4])),1]
      countm1<-apply(input[,grep(".6P|X6.P",colnames(input))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
      countm2<-apply(input[,grep(".7P|X7.P",colnames(input))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
      countm3<-apply(input[,grep("NC.P",colnames(input))[Samping[[i]]$test]],2,function(x) tapply(x,bio$V5,function(x) sum(x>best,na.rm=T)))
      Lnum1<-c(Lnum1,c=sum(apply(countm1,2,function(x) which.max(x)==2)))
      Lnum2<-c(Lnum2,c=sum(apply(countm2,2,function(x) which.max(x)==6)))
      Lnum3<-c(Lnum3,c=sum(apply(countm3,2,function(x) which.max(x)==10)))
      
      acc1<-rbind(acc1,c(best,sum(apply(countm1,2,function(x) which.max(x)==2))/(length(Samping[[i]]$test))))
      acc2<-rbind(acc2,c(best,sum(apply(countm2,2,function(x) which.max(x)==6))/(length(Samping[[i]]$test))))
      acc3<-rbind(acc3,c(best,sum(apply(countm3,2,function(x) which.max(x)==10))/(length(Samping[[i]]$test))))
      
      Best<-c(Best,best)
    }
    cc1<-apply(test[,grep(".6P|X6.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
    cc2<-apply(test[,grep(".7P|X7.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
    cc3<-apply(test[,grep("NC.P",colnames(test))],2,function(x) tapply(x,bio$V5,function(x) sum(x>mean(Best),na.rm=T)))
    ccc1=sum(apply(cc1,2,function(x) which.max(x)==2))/5
    ccc2=sum(apply(cc2,2,function(x) which.max(x)==6))/5
    ccc3=sum(apply(cc3,2,function(x) which.max(x)==10))/5
    ccc1
    ccc2
    ccc3
    rlt<-c(mean(Best),ccc1,ccc2,ccc3)
    Rlt<-rbind(Rlt,rlt)
  }
  return(Rlt)
}
