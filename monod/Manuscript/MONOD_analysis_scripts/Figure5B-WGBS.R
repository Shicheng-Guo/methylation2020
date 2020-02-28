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

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
}

gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  GSI<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-c()
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    if(length(gmax)<1){print(data[i,])}
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
      gsit<-c(gsit,tmp)
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    GSI<-c(GSI,sum(gsit,na.rm=T))
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

topgsi2bio<-function(topgsi){
  cor2bed<-function(cor){
    cor<-as.character(cor)
    a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
    bed<-matrix(a,ncol=3,byrow=T)
    return(data.frame(bed))
  }
  bio<-data.frame(cor2bed(topgsi[,1]),topgsi[,2:3])
  rownames(bio)<-topgsi[,1]
  return(bio)
}

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}   

rename<-function(data){
  Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P",colnames(data))]
  colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
  colnames(Data)<-gsub("[.]","-",colnames(Data))
  colnames(Data)[grep("age|new|centenarian",colnames(Data))]<-"WBC"
  colnames(Data)[grep("X6.T|CTT|CCT",colnames(Data))]<-"CCT"
  colnames(Data)[grep("X7.T|LCT",colnames(Data))]<-"LCT"
  colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
  colnames(Data)[grep("X6.P|RRBS.6P",colnames(Data))]<-"CCP"
  colnames(Data)[grep("X7.P|RRBS.7P",colnames(Data))]<-"LCP"
  colnames(Data)[grep("NC.P",colnames(Data))]<-"NP"
  colnames(Data)[grep("X7.T|SRX381722|SRX381719|SRX381716",colnames(Data))]<-"LCT"
  colnames(Data)[grep("X6.T|CTT|SRX381569",colnames(Data))]<-"CCT"
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

OneTimeTissuePrediction<-function(data,bio){
  input<-data[match(rownames(bio),rownames(data)),]
  Num<-c()
  for(j in seq(0,0.6,0.05)){
    counts1<-apply(input[,grep(".6P|X6.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts2<-apply(input[,grep(".7P|X7.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts3<-apply(input[,grep("NC.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    num<-data.frame(id=j,c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                    c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                    c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  print(counts1)
  acc<-sweep(Num, 2,c(1,30,29,75) , `/`)
  rlt<-acc[which.max(rowSums(acc[,2:4])),]
  return(rlt)
}


GsirltfilterBinaryMatrix<-function(DATA,Gsirlt){
  
  input<-DATA[match(Gsirlt[,1],rownames(DATA)),]
  name<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[.]"))[1]))
  colnames(input)<-name
  Mean<-data.frame(apply(input,1,function(x) tapply(x,name, function(x) mean(x,na.rm=T))))
  
}

Gsirltfilter<-function(DATA,Gsirlt,thres){
  inc<-c()
  for(i in 1:nrow(Gsirlt)){
    loc<-Gsirlt[i,1]
    group<-Gsirlt[i,2]
    input<-DATA[match(loc,rownames(DATA)),grep(group,colnames(DATA))]
    Mean1<-mean(as.numeric(input),na.rm=T)
    xinput<-na.omit(unlist(DATA[match(loc,rownames(DATA)),(colnames(DATA)!=group)]))
    names(xinput)<-unlist(lapply(names(xinput),function(x) unlist(strsplit(x,"[.]"))[1]))
    Mean2<-tapply(as.numeric(xinput),names(xinput),function(x) mean(x,na.rm=T))
    if(Mean1>thres && all(Mean2<thres)){
      inc<-c(inc,i)
    }
    print(i)
  }
  Gsirlt<-Gsirlt[inc,]
  return(Gsirlt)
}

reorderdata<-function(data,referenceorder){
  # reference<-c("Brain","Colon","Intestine","Kidney","Liver","Lung","Pancreas","Spleen","Stomach","WBC","CCT","CCP","LCT","LCP","NP")
  column<-c()
  for(i in referenceorder){
    column<-c(column,grep(i,colnames(data)))
  }
  output<-data[,column]
  colnames(output)<-unlist(lapply(colnames(output),function(x) unlist(strsplit(x,"[.]"))[1]))
  return(output)
}

##################################################################################
####################### DATA LOAD AND PROCESS ################################
##################################################################################
setwd("/home/shg047/oasis/monod/hapinfo")
saminfo<-read.table("/home/shg047/oasis/monod/saminfo/N37Salk.saminfo",sep="\t")
data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
# load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
# feature reduction (WGBS+RRBS missing<60, each plasma category(colon,lung and normal)<60%, Plasma missing<50%)
rm<-frm(data)
data<-data[-rm,]
#input<-data[match(rownames(bio),rownames(data)),]
# copy our biomarker form supplementary table to: /home/sguo/Dropbox/Project/methylation/monod/Manuscript/MONOD_analysis_scripts/biomarker2.txt (download from supplementary)
#bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)<-gsub("[.]","-",colnames(Data))
colnames(Data)[grep("age|new|centenarian|WB",colnames(Data))]<-"WBC"
colnames(Data)[grep("X7.T|X6.T|SRX|CTT",colnames(Data))]<-"CT"
colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
# for tissue-specific biomarkers
library("impute")
library("preprocessCore")
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
#DATA<-Data[match(rownames(bio),rownames(Data)),grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
colnames(DATA)<-colnames(Data)[grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))]
#factor<-max(apply(DATA,2,function(x) mean(x,na.rm=T)))/apply(DATA,2,function(x) mean(x,na.rm=T))
#DATA<-sweep(DATA, 2,factor, `*`)
#for(i in names(table(colnames(DATA)))){
#  DATA[,colnames(DATA)==i]<-normalize.quantiles(data.matrix(DATA[,colnames(DATA)==i]))
#}
#DATA<-normalize.quantiles(data.matrix(DATA))
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
DATA<-RawNARemove(DATA,missratio=0.60)
#DATA1<-ColNARemove(DATA,missratio=0.79)
#DATA2<-impute.knn(data.matrix(DATA1))$data
#gsirlt<-gsi(DATA2)
##################################################################################
####################### GSI AND FURTHER SELECTION ################################
##################################################################################
colnames(DATA)
gsirlt<-gsi(DATA)
save.image("Figure5A-WGBS-biomarker.RData")
load("Figure5A-WGBS-biomarker.RData")

# check the heatamp for high gsi
setwd("/home/shg047/oasis/monod/hapinfo")
bio<-topgsi2bio(TopGSIByCategory(gsi(DATA),top=300))
saminfo<-read.table("/home/shg047/oasis/monod/saminfo/N37Salk.saminfo",sep="\t")
data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T|NC.P|7.P|RRBS.7P|6.P|RRBS.6P|CCT|LCT",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)<-gsub("[.]","-",colnames(Data))
colnames(Data)[grep("age|new|centenarian|WB|WBC",colnames(Data))]<-"WBC"
colnames(Data)[grep("X7.T|SRX381722|SRX381719|SRX381716",colnames(Data))]<-"LCT"
colnames(Data)[grep("X6.T|CTT|SRX381569",colnames(Data))]<-"CCT"
colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
colnames(Data)
input<-data[match(rownames(bio),rownames(data)),]
input<-rename(input)
colnames(input)
reference<-c("Brain","Colon","Intestine","Kidney","Liver","Lung","Pancreas","Spleen","Stomach","WBC","CCT","CCP","LCT","LCP","NP")
input<-reorderdata(input,reference)
colnames(input)

heatmap<-function(mydata,output){
  tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  col1<-tol21rainbow[as.numeric(as.factor(rownames(mydata)))]
  col2<-tol21rainbow[as.numeric(as.factor(colnames(mydata)))]
  col=colorRampPalette(c("yellow", "blue"))(20) 
  require(gplots)    
  outputfile<-paste(output,".pdf",sep="")
  pdf(outputfile)
  par(mar=c(5,1,5,0))
  heatmaprlt<-heatmap.2(as.matrix(mydata),hclustfun=hclustfunc, distfun=distfunc,
                        RowSideColors=col1, 
                        ColSideColors=col2,
                        labRow=F,
                        trace="none",
                        Colv=F,
                        Rowv=F,
                        col=col,
                        density.info="none")
  legend=unique(data.frame(col2,as.factor(colnames(mydata))))
  legend(x=0.9,y=0.8,legend=legend[,2],col=as.character(legend[,1]),pch=15,cex = 0.5)
  #legend=unique(data.frame(col1,as.factor(rownames(mydata))))
  #legend(x=0.9,y=0.5,legend=legend[,2],col=as.character(legend[,1]),pch=15,cex = 0.5)
  dev.off()
}
heatmap(input,output="Figure-HighGSI-in-solid-tissue-and-plamsma.pdf")

# extent feature selection(variable selection, data driven)

NoneWBCData<-DATA[match(subset(gsirlt,group=="WBC"|group=="Spleen")[,1],rownames(DATA)), colnames(DATA) !="WBC" & colnames(DATA) !="Spleen"]
Gsirlt<-gsirlt[-match(rownames(NoneWBCData[which(apply(NoneWBCData,1,function(x) any(x>0.004,na.rm=T))),]),gsirlt$region),]
#bio1<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
DATA2<-DATA[na.omit(match((subset(Gsirlt,group!="WBC" & group!="Spleen")[,1]),rownames(DATA))),grep("WBC|Spleen",colnames(DATA))]
DATA3<-DATA[na.omit(match((subset(Gsirlt,group="WBC"|group!="Spleen")[,1]),rownames(DATA))),colnames(DATA)!="WBC" &colnames(DATA)!="Spleen"]
Acc<-c()
for(i in seq(0.01,1,0.01)){
rowmean1<-names(which(unlist(apply(DATA2,1,function(x) any(x>0.1)))))
rowmean2<-names(which(unlist(apply(DATA3,1,function(x) any(x>0.001)))))
DATA4<-DATA[-match(unique(c(rowmean1,rowmean2)),rownames(DATA)),]
Gsirlt<-gsirlt[unique(na.omit(match(rownames(DATA4),gsirlt[,1]))),]
#Gsirlt<-Gsirltfilter(DATA4,Gsirlt,thres=0.05)
topgsi<-TopGSIByCategory(Gsirlt,top=500)   
bio<-topgsi2bio(topgsi)
acc<-OneTimeTissuePrediction(data,bio)
acc
Acc<-rbind(Acc,acc)
print(round(unlist(c(i,acc)),2))
}

# take all the biomarkers
Gsirlt<-Gsirlt[na.omit(match(rownames(DATA),Gsirlt[,1])),]
xx1<-gsirlt[gsirlt[,1] %in% rownames(bio1),]
xx2<-bio1[na.omit(match(gsirlt[,1],rownames(bio1))),]
sum(is.na(match(rownames(bio),Gsirlt[,1])))
# compared with our previsou identified biomarker
Bio<-bio[sample(1:nrow(bio),0.9*nrow(bio)),]  # apply parts of the tissue-specific biomarkers so that we can evalate the stability
Bio<-Bio[order(Bio[,4]),]
DATA3<-DATA[match(rownames(Bio),rownames(DATA)),grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WB",colnames(DATA))] 
matchrlt<-apply(DATA3,2,function(x) names(table(Bio[,4]))[which.max(tapply(x,Bio[,4],function(x) sum(x>0.1,na.rm=T)))])   # only need to set a small threshold to control the background noise since it is not plasma data.
for(i in names(table(Bio[,4]))){
  acc=sum(matchrlt[grep(i,names(matchrlt))]==i)/length(grep(i,names(matchrlt)))
  print(c(i,acc))
}

print(c(i,length(names(which(!apply(MM,1,function(x) any(x>0.004,na.rm=T)))))))
apply(DATA[which(!apply(MM,1,function(x) any(x>0.004,na.rm=T))), colnames(DATA) =="WBC"],2,function(x) range(x,na.rm=T))
sum(names(which(!apply(MM,1,function(x) any(x>0.0004,na.rm=T)))) %in% rownames(subset(bio,V5=="WBC")))
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
bio[,6]=rownames(bio)
x=merge(bio,gsirlt,by.x=6,by.y=1)
sum(x$V5==x$group)
y=(subset(x,V5=="Colon"))
table(y$group)
y=(subset(x,V5=="WBC"))
DATA[match(y[,1],rownames(DATA)),grep("WBC",colnames(DATA))]
table(y$group)
y=(subset(x,V5=="Lung"))
table(y$group)


######################################################################
##################### One-time test ##################################
######################################################################
# The most easy situation is how the performance for one-time prediction
# Take different `t` and then use all the samples, give the prediction arracy
# Response: Random select 5 samples for each group, leave them alone and then apply our modle to these samp
pdf("one-time-prediction-test-acc.pdf")
plot(c1~id,acc,type="l",col=2,cex=2,lwd=2,ylim=c(0,1),xlab="MHL threshold (t)",ylab="Accuracy of the prediction")
lines(c2~id,acc,lty=2,col=3,cex=2,lwd=2)
lines(c3~id,acc,lty=3,col=4,cex=2,lwd=2)
legend("topright",legend=c("CCP","LCP","NCP"),lty=c(1,2,3),col=c(2,3,4),lwd=2, bty="n")
dev.off()
print(ncc)


