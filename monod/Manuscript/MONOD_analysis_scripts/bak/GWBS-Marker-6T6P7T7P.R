# 2017-05-10
# Tissue specific MHL markers could predict 6T, 6P, 7T, 7P with 100% accuracy.

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("impute")  
#install.packages("gplots")
#install.packages("RColorBrewer")
#install.packages("grDevices")

library("gplots")
library("RColorBrewer")
library("grDevices")
library("impute")
library("reshape2")
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
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    if(length(gmax)<1){print(data[i,])}
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
      gsit<-c(gsit,tmp)
    }
    print(i)
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
  Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|WB|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P",colnames(data))]
  colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
  colnames(Data)<-gsub("[.]","-",colnames(Data))
  colnames(Data)[grep("age|new|centenarian|middle",colnames(Data))]<-"WBC"
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

rbedintersect<-function(bed1,ref){
  Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
    #create temp files
    a.file=tempfile()
    b.file=tempfile()
    out   =tempfile()
    options(scipen =99) # not to use scientific notation when writing out
    #write bed formatted dataframes to tempfile
    write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
    # create the command string and call the command using system()
    command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
    cat(command,"\n")
    try(system(command))
    res=read.table(out,header=F)
    unlink(a.file);unlink(b.file);unlink(out)
    res=subset(res,V5!=".")
    return(res)
  }
  merge<-Rbedtools(functionstring="intersectBed",bed1,ref,opt.string="-wao")
  return(merge)
}

bed2cg<-function(bed1){
  ref<-read.table("~/work/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
  cor2bed<-function(cor){
    cor<-as.character(cor)
    a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
    bed<-matrix(a,ncol=3,byrow=T)
    bed<-data.frame(bed,cor)
    return(data.frame(bed))
  }
  rbedintersect<-function(bed1,ref){
    Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
      #create temp files
      a.file=tempfile()
      b.file=tempfile()
      out   =tempfile()
      options(scipen =99) # not to use scientific notation when writing out
      #write bed formatted dataframes to tempfile
      write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
      write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
      # create the command string and call the command using system()
      command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
      cat(command,"\n")
      try(system(command))
      res=read.table(out,header=F)
      unlink(a.file);unlink(b.file);unlink(out)
      res=subset(res,V5!=".")
      return(res)
    }
    merge<-Rbedtools(functionstring="intersectBed",bed1,ref,opt.string="-wao")
    return(merge)
  }
  merge<-rbedintersect(bed1,ref)
  return(merge)
}
cor2bed<-function(cor){
  cor<-as.character(cor)
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  bed<-data.frame(bed,cor)
  return(data.frame(bed))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
OneTimeTissuePredictionTrainStage<-function(data,bio){
  input<-data[match(rownames(bio),rownames(data)),]
  Num<-c()
  for(j in seq(0,0.9,0.05)){
    counts1<-apply(input[,grep(".6P|X6.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts2<-apply(input[,grep(".7P|X7.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts3<-apply(input[,grep("NC.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    num<-data.frame(id=j,c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                    c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                    c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  print(counts1)
  acc<-sweep(Num, 2,c(1,5,5,5) , `/`)
  rlt<-acc[which.max(rowSums(acc[,2:4])),]
  return(rlt)
}
OneTimeTissuePredictionTrainStage<-function(data,bio,freq){
  input<-data[match(rownames(bio),rownames(data)),]
  Num<-c()
  for(j in seq(0,0.9,0.05)){
    counts1<-apply(input[,grep(".6P|X6.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts2<-apply(input[,grep(".7P|X7.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts3<-apply(input[,grep("NC.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts1<-sweep(counts1, 1,table(bio[,4]), `/`)
    counts2<-sweep(counts2, 1,table(bio[,4]), `/`)
    counts3<-sweep(counts3, 1,table(bio[,4]), `/`)
    num<-data.frame(id=j,c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                    c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                    c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  print(counts1)
  acc<-sweep(Num, 2,c(1,freq) , `/`)
  rlt<-acc[which.max(rowSums(acc[,2:4])),]
  return(rlt)
}
cg2bed<-function(cg,extend=0){
  bed2cor<-function(bed){
    cor<-apply(bed,1,function(x) paste(x[1],":",as.numeric(x[2])-extend,"-",as.numeric(x[3])+extend,sep=""))
    cor<-gsub(" ","",cor)
    return(cor)
  }
  ref<-read.table("~/work/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
  bed<-ref[match(cg,ref[,4]),1:3]
  bed[,2]=bed[,2]-extend
  bed[,3]=bed[,3]+extend
  cor<-bed2cor(bed)
  rlt<-data.frame(bed,cor,cg)
  return(rlt)
}



##################################################################################
####################### DATA LOAD AND PROCESS ################################
##################################################################################
setwd("/home/shg047/oasis/monod/hapinfo")
saminfo<-read.table("/home/shg047/oasis/monod/saminfo/N37Salk.saminfo",sep="\t")
#bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
#data<-read.table("MHL4.txt",head=T,row.names=1,sep="\t")
#save(data,file="MHL4.RData")
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
# feature reduction (WGBS+RRBS missing<60, each plasma category(colon,lung and normal)<60%, Plasma missing<50%)
#rm<-frm(data)
#data<-data[-rm,]
Data=data[,grep("STL|N37|age|ENC|new|centenarian|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
Data<-rename(Data)
colnames(Data)<-unlist(lapply(colnames(Data),function(x) unlist(strsplit(x,"[.|-]"))[1]))
#DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC|CCT|CCP|LCT|LCP|NP",colnames(Data))] 
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
gsirlt<-gsi(DATA)

colon<-subset(gsirlt,group=="Colon")
head(colon[order(colon[,3],decreasing=T),])
data[match("chr7:127992139-127992150",rownames(data)),]
DATA[match("chr7:127992139-127992150",rownames(DATA)),]
Data[match("chr7:127992139-127992150",rownames(Data)),]
save(gsirlt,file="gsirlt.RData")
system("scp /oasis/tscc/scratch/shg047/monod/hapinfo/gsirlt.RData shg047@genome-miner.ucsd.edu:/media/NAS1/shg047/tcga/pancancer")

#### Go to genome-miner ######
setwd("/media/NAS1/shg047/tcga/pancancer")
load("gsirlt.RData")
load("PancancerMethMatrix_March2016.RData")
colon<-subset(gsirlt,group=="Colon")
# make sure colon specfici markers hypemethyation
cg<-bed2cg(cor2bed(colon[,1]))
input<-data[unique(match(cg[,8],rownames(data))),grep("COAD",colnames(data))]
input<-input[which(apply(input,1,function(x) mean(x,na.rm=T))>0.2),]
bed<-unique(rbedintersect(cor2bed(colon[,1]),unique(cg2bed(rownames(input))[,1:4]))[,1:4])
system("scp /media/NAS1/shg047/tcga/pancancer/colon.450h.RData shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/monod/hapinfo")
system("scp /media/NAS1/shg047/tcga/pancancer/PancancerMethMatrix_March2016.RData shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/monod/hapinfo")
# make sure these biomarker is low methylated in WGBS
setwd("Downloads")
load("colon.450h.RData")
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC|CCP|LCP|NP",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
test<-DATA[match(bed[,4],rownames(DATA)),]
wbcexclude1<-which(apply(test[,grep("WBC",colnames(test))],1,function(x) mean(x,na.rm=T))>0.001)
wbcexclude2<-which(apply(test[,grep("Spleen",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude3<-which(apply(test[,grep("Stomach",colnames(test))],1,function(x) mean(x,na.rm=T))>0.001)
wbcexclude4<-which(apply(test[,grep("Kidney",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude5<-which(apply(test[,grep("Intestine",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude6<-which(apply(test[,grep("Liver",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude7<-which(apply(test[,grep("Pancreas",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude8<-which(apply(test[,grep("Lung",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude9<-which(apply(test[,grep("Brain",colnames(test))],1,function(x) mean(x,na.rm=T))>0.01)
wbcexclude<-unique(c(wbcexclude1,wbcexclude2,wbcexclude3,wbcexclude4,wbcexclude5,wbcexclude6,wbcexclude7,wbcexclude8,wbcexclude9))
wbcexclude<-unique(c(wbcexclude1))
test<-test[-wbcexclude,]
xx<-apply(test,2,function(x) sum(na.omit(x)>0.005)/length(na.omit(x)))
xx<-data.frame(mdf=xx,type=names(xx))
xx$type <- factor(xx$type,levels = c("NP","LCP","CCP","WBC","Brain","Colon","Intestine",
                                     "Kidney","Liver","Lung","Pancreas","Spleen","Stomach"
                                     ),ordered = TRUE)
uxx<-melt(test)
colnames(uxx)<-c("type","mhl")
uxx$type <- factor(uxx$type,levels = c("NP","LCP","CCP","WBC","Brain","Colon","Intestine",
                                     "Kidney","Liver","Lung","Pancreas","Spleen","Stomach"
                                     ),ordered = TRUE)
png("assign.number.png")
p1<-ggplot(aes(y = mhl, x = type, fill = type), data = uxx) + 
  geom_boxplot(outlier.shape =NA,outlier.colour="blue")+ 
  coord_flip()+
  geom_point(position = position_jitter(width = 0.2),size=0.5)

p2<-ggplot(aes(y = mdf, x = type, fill = type), data = xx) + 
  geom_boxplot(outlier.shape =NA,outlier.colour="blue")+ 
  coord_flip()+
  geom_point(position = position_jitter(width = 0.2),size=0.5)
multiplot(p1, p2, cols=2)
dev.off()





library("ggplot2")
png("assign.number.png")
#boxplot(xx~names(xx),)
ggplot(aes(y = mhl, x = type, fill = type), data = xx) + geom_boxplot(outlier.shape =NA,outlier.colour="blue")+ coord_flip()
dev.off()
# 
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC|CCP|LCP|NP",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
xx<-apply(DATA,1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T)))
yy<-data.frame(t(xx))
zz<-subset(yy,Colon>0.3 & WBC<0.03 & Lung<0.03)
dim(zz)
test<-DATA[match(rownames(zz),rownames(DATA)),]
xx<-apply(test,2,function(x) sum(na.omit(x)>0.01)/length(na.omit(x)))
xx
png("assign.number.png")
boxplot(xx~names(xx))
dev.off()


# check GSI which are not match and correct them with our newdata
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
DATA<-DATA[match(rownames(bio),rownames(DATA)),]
gsirlt<-gsi(DATA)
bio2<-data.frame(bio,gsirlt[match(rownames(bio),gsirlt[,1]),])
table(bio2$V5,bio2$group)

xx1<-bio2[which(bio2$V5=="Colon" & bio2$group!="Colon"),]
yy1<-apply(DATA,1,function(x) tapply(x,colnames(DATA), function(x) mean(x,na.rm=T)))
yy1<-yy1[,match(rownames(xx1),colnames(yy1))]
yy1[,1:5]
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC|CCP",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
yy1<-apply(DATA,1,function(x) tapply(x,colnames(DATA), function(x) mean(x,na.rm=T)))
yy1<-yy1[,match(rownames(xx1),colnames(yy1))]
yy1[,1:5]


xx2<-bio2[which(bio2$V5=="Lung" & bio2$group!="Lung"),]
yy2<-apply(DATA,1,function(x) tapply(x,colnames(DATA), function(x) mean(x,na.rm=T)))
yy2<-yy2[,match(rownames(xx1),colnames(yy1))]
yy2[,1:5]
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC|LCP",colnames(Data))] 
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
yy2<-apply(DATA,1,function(x) tapply(x,colnames(DATA), function(x) mean(x,na.rm=T)))
yy2<-yy2[,match(rownames(xx2),colnames(yy2))]
yy2[,1:5]

match(rownames(xx2),colnames(apply(DATA,1,function(x) tapply(x,colnames(DATA), function(x) mean(x,na.rm=T)))))
match(rownames(xx3),colnames(apply(DATA,1,function(x) tapply(x,colnames(DATA), function(x) mean(x,na.rm=T)))))
match(rownames(xx4),colnames(apply(DATA,1,function(x) tapply(x,colnames(DATA), function(x) mean(x,na.rm=T)))))

gsi1<-function(DATA,bio){
  rlt<-list()
  grep("LCT",colnames(DATA))
  DATA<-DATA[,-grep("CCP|CCT",colnames(DATA))]
  colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
  rlt$data=DATA
  rlt$element<-apply(DATA[which(bio[,4]=="Colon"),],1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T)))
  rlt$rlt<-apply(rlt$element,2,function(x) which.max(x))
  rlt$xtab<-table(rlt$rlt)
  rlt$C2L<-DATA[match(names(which(apply(apply(DATA[which(bio[,4]=="Colon"),],1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T))),2,function(x) which.max(x))==8)),rownames(DATA)),grep("LCT",colnames(DATA))]
  rlt$mmax<-apply(rlt$C2L,2,function(x) mean(x,na.rm=T))
  rlt$whichmax<-apply(DATA,1,function(x) order(x,decreasing=T))
  return(rlt)
}

gsi2<-function(DATA,bio){
rlt<-list()
grep("LCT",colnames(DATA))
#DATA<-DATA[,-grep("CCP",colnames(DATA))]
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
rlt$element<-apply(DATA[which(bio[,4]=="Colon"),],1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T)))
rlt$rlt<-apply(rlt$element,2,function(x) which.max(x))
rlt$xtab<-table(rlt$rlt)
rlt$C2L<-DATA[match(names(which(apply(apply(DATA[which(bio[,4]=="Colon"),],1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T))),2,function(x) which.max(x))==8)),rownames(DATA)),grep("LCT",colnames(DATA))]
rlt$mmax<-apply(rlt$C2L,2,function(x) mean(x,na.rm=T))
return(rlt)
}
test1=data[,c(grep("X7.T|X6.T",colnames(data)),grep("NC.P",colnames(data))[1:10],grep("X6.P",colnames(data))[c(1,3,4,5,6)],grep("X7.P",colnames(data))[c(1,3,4,5,6)])]
test2=data[,c(grep("NC.P",colnames(data))[1:5],grep("X6.P",colnames(data))[c(1,3,4,5,6)],grep("X7.P",colnames(data))[c(1,3,4,5,6)])]
test3=data[,c(grep("NC.P",colnames(data))[1:75],grep("X6.P|RRBS.6P",colnames(data))[c(1:30)],grep("X7.P|RRBS.7P",colnames(data))[c(1:29)])]
OneTimeTissuePredictionTrainStage(test1,bio2,c(10,10,10))
OneTimeTissuePredictionTrainStage(test2,bio2,c(5,5,5))
OneTimeTissuePredictionTrainStage(test3,bio2,c(30,29,75))



table(apply(ccrlt1$whichmax,2,function(x) colnames(ccrlt1$data)[x[1]]))
table(apply(ccrlt1$whichmax,2,function(x) colnames(ccrlt1$data)[x[2]]))

ccrlt1<-gsi1(DATA,bio)
ccrlt1$xtab
ccrlt2<-gsi2(DATA,bio)
ccrlt2$xtab

write.table(cor2bed(names(which(ccrlt2$rlt==2))),file="Colon77CCP.txt",sep="\t",col.names=F,row.names=F,quote=F)
DATA[match(names(which(ccrlt2$rlt==2)),rownames(DATA)),grep("Colon",colnames(DATA))]
table(apply(apply(DATA[which(bio[,4]=="Lung"),],1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T))),2,function(x) which.max(x)))
C2L<-DATA[match(names(which(apply(apply(DATA[which(bio[,4]=="Lung"),],1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T))),2,function(x) which.max(x))==3)),rownames(DATA)),grep("LCT",colnames(DATA))]
apply(C2L,2,function(x) mean(x,na.rm=T))

bio<-read.table("~/work/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
rlt<-bed2cg(cor2bed(rownames(subset(bio,V5=="Colon"))))
xdata<-data[match(rlt$V8,rownames(data)),grep("COAD",colnames(data))]
length(which(apply(xdata,1,function(x) mean(x,na.rm=T))>0.1))

rlt2<-bed2cg(read.table("Colon77CCP.txt"))
xdata2<-data[match(rlt2$V8,rownames(data)),grep("COAD",colnames(data))]
length(which(apply(xdata2,1,function(x) mean(x,na.rm=T))>0.1))



