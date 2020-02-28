---
title: "MONOD Project Discovery"
author: "Shicheng Guo"
date: "Wednesday, March 30, 2016"
output: html_document
---


```{r}
data<-read.table("monod.mhl.List1.march29.txt",head=T,row.names=1,sep="\t",as.is=T,check.names=F)
sum(is.na(data))/(nrow(data)*ncol(data))

# re-group/collect the samples #colon cancer
Group<-colnames(data)
Group[grep("CTR|NC|Pregn",Group)]<-"NP"
Group[grep("6-P-",Group)]<-"CCP"
Group[grep("6-T-|CTT-|metastasis_colon|tumor_colon",Group)]<-"CCT"
Group[grep("normal_colon|N37-Colon|SG-01",Group)]<-"NCT"
Group[grep("WB-",Group)]<-"WB"
Group[grep("STL|N37-|methylC-",Group)]<-"ONT"
newdata<-data
colnames(newdata)<-Group
newdata<-data.frame(newdata,check.names=F)
YY<-grep("NP|CCP|CCT|NCT|WB|ONT",Group)
newdata<-newdata[,YY]
idx<-sapply(colnames(newdata),function(x) unlist(strsplit(x,"[.]"))[1])
colnames(newdata)<-idx
newdata<-data.matrix(newdata)

# feature selection # 1: CT> CP/NT and NT,NP,ONT<0.1
rlt<-apply(newdata,1,function(x) tapply(x,idx,function(x) median(x,na.rm=T)))
target<-c()
for(i in 1:ncol(rlt)){
  order=order(rlt[,i],decreasing=T)
  if(order==c(2,1,3,4,5,6) ||order==c(2,3,1,4,5,6)){
    if(! all(is.na(rlt[1,i])) & ! all(is.na(rlt[4,i])) & ! all(is.na(rlt[5,i])) & ! all(is.na(rlt[6,i])) & rlt[1,i]>0.1 & rlt[4,i]<0.1 & rlt[5,i]<0.1 & rlt[6,i]<0.1 ){
      print (i)
      target<-c(target,i)
    }
  } 
}


# feature selection # 2: CT>0.2 and NP<0.05 and WB<0.05
# rlt<-apply(newdata,1,function(x) tapply(x,idx,function(x) mean(x,na.rm=T)))
target<-c()
for(i in 1:ncol(rlt)){
  if(! all(is.na(rlt[2,i])) & ! all(is.na(rlt[4,i])) & ! all(is.na(rlt[6,i])) & rlt[2,i]>0.2 & rlt[4,i]<0.1 & rlt[6,i]<0.1){
    print (i)
    target<-c(target,i)
  } 
}

# feature selection # 2: CT>0.5 and WB<0.1
# rlt<-apply(newdata,1,function(x) tapply(x,idx,function(x) mean(x,na.rm=T)))
target<-c()
for(i in 1:ncol(rlt)){
  if(! all(is.na(rlt[3,i])) & ! all(is.na(rlt[4,i])) & ! all(is.na(rlt[5,i])) & ! all(is.na(rlt[6,i])) & rlt[3,i]>0.5 & rlt[4,i]<0.1 & rlt[5,i]<0.1 & rlt[6,i]<0.1){
    print (i)
    target<-c(target,i)
  } 
}
length(target)
# ggplot for boxplot
newdataGGplot<-newdata[match(colnames(rlt)[target],rownames(newdata)),]
write.table(newdataGGplot,file="monond-colon-plasma.txt",col.names=NA,row.names=T,sep="\t",quote=F)
```

You can also embed plots, for example:

```{r, echo=FALSE}
plot(cars)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
