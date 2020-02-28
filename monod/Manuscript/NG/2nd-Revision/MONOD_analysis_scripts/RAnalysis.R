# This software is Copyright ? 2017 The Regents of the University of California. All Rights Reserved.
#  
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
#  
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# kzhang@ucsd.edu

##===============================================================================================================================================================
##   Summarizes data and analysis for MONOD project
##   Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples and tissue-of-origin mapping from plasma DNA
##   Functions for MHB, MHL, GSI, heatmap, deconvolution, tissue-of-origin mapping, 
##================================================================================================================================================================

#!/usr/bin/R
# Methylation haplotype and Methylation haplotype load (mode weigth i)
# Version 1.3
# Update: Jan/9/2017


## Figure 1
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
  return(res)
}
Rbedtoolsort<-function(functionstring="bedtools sort",bed1,opt.string=""){
  #create temp files
  a.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-i",a.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}
bedwithgap<-function(bed,gap){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
}
cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

bed2cor<-function(bed){
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],"-",unlist(strsplit(x,"\t"))[3],sep="")})
  return(cor)
}
bedEnrichmentTest<-function(bed1,bed2,assembly="hg19",iteration=1000000){
  rlt<-list()
  hsa.chr<-paste("chr",c(1:22,"X","Y"),sep="")
  hg19.size<-c("249250621","243199373","198022430","191154276","180915260","171115067",
               "159138663","146364022","141213431","135534747","135006516","133851895",
               "115169878","107349540","102531392","90354753","81195210","78077248",
               "59128983","63025520","48129895","51304566","155270560","59373566")
  hg38.size<-c("248956422","242193529","198295559","190214555","181538259","170805979",
               "159345973","145138636","138394717","133797422","135086622","133275309",
               "114364328","107043718","101991189","90338345","83257441","80373285",
               "58617616","64444167","46709983","50818468","156040895","57227415")
  mm.chr<-paste("chr",c(1:19,"X","Y"),sep="")
  mm9.size<-c("197195432","181748087","159599783","155630120","152537259","149517037",
              "152524553","131738871","124076172","129993255","121843856","121257530",
              "120284312","125194864","103494974","98319150","95272651","90772031",
              "61342430","166650296","15902555")
  mm10.size<-c("195471971","182113224","160039680","156508116","151834684","149736546",
               "145441459","129401213","124595110","130694993","122082543","120129022","120421639",
               "124902244","104043685","98207768","94987271","90702639","61431566","171031299","91744698")
  hg19.chrom.sizes<-data.frame(hsa.chr,hg19.size)
  hg38.chrom.sizes<-data.frame(hsa.chr,hg38.size)
  mm9.chrom.sizes<-data.frame(mm.chr,mm9.size)
  mm10.chrom.sizes<-data.frame(mm.chr,mm10.size)
  write.table(hg19.chrom.sizes,file="hg19.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(hg38.chrom.sizes,file="hg38.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(mm9.chrom.sizes,file="mm9.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(mm10.chrom.sizes,file="mm10.chrom.sizes",col.names=F,row.names=F,quote=F,sep="\t")
  
  cor2bed<-function(cor){
    a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
    bed<-matrix(a,ncol=3,byrow=T)
    return(data.frame(bed))
  }
  
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
    return(res)
  }
  
  command=paste("bedtools intersect -wa -u","-a",bed1,"-b",bed2,">","output1",sep=" ")
  try(system(command))
  count1<-read.table("output1")
  observe<-nrow(count1)
  try(system(paste("rm","output1",sep=" " )))
  try(system(paste("rm","output2",sep=" " )))
  chrom.size=paste(assembly,"chrom.sizes",sep=".")
  for (i in 1:iteration){
    command=paste("bedtools shuffle","-i",bed1,"-g",chrom.size,"| bedtools intersect -wa -u -a -","-b",bed2,"|wc -l >> output2",sep=" ")
    try(system(command))
    print(i)
  }
  count2<-read.table("output2")
  
  print(paste("The expected overlap region number is",mean(count2[,1]),sep=" "))
  print(paste("The observed overlap region number is",observe,sep=" "))
  
  pvalue<-1-sum(count2[,1] < observe)/iteration
  if(pvalue==0){
    pvalue=paste("<",1/iteration,sep=" ");
  }else{
    pvalue=paste("=",pvalue,sep=" ")
  }
  print(paste("The P-value of the target bed region overlapped with source bed, Pvalue",pvalue,sep=""))
  try(system("mv output2 expected.overlab.counts.number"))
  
  pdf("enrichment.hist.box.pdf")
  hist(count2[,1],breaks=25,col="red",main="",
       ylab="Counts")
  boxplot(count2[,1], at=15,horizontal=TRUE,
          outline=F,boxwex=50, frame=F,border=3,col = "green1", add = TRUE,lwd=6)
  dev.off()
  
  rlt$expectation<-count2[,1]
  rlt$observation<-observe
  rlt$pvalue=pvalue
  return(rlt)
}

## Figure 2
hap_comb_tmp<-function(hapinfo.single){
  # hap_com.single is a character of haplotype, such as "TCTCTT"
  hap_com.single<-c()
  for(i in 1:nchar(hapinfo.single)){
    for(j in i:(nchar(hapinfo.single))){
      hap_com.single<-c(hap_com.single,(substr(hapinfo.single,i,j)))
    } 
  }
  return(hap_com.single)
}

hap_comb<-function(hapinfo){
  # haplotype is array of observed methylation haplotype and return
  hap_comb<-unlist(lapply(hapinfo,hap_comb_tmp))
  return(hap_comb)
}

mhl<-function(hap_comb_array){
  mhl_value<-c()
  hap_comb_array<-unlist(hap_comb_array)
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  unmeth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="T"
  hap_comb_array<-c(hap_comb_array[meth_hap],hap_comb_array[unmeth_hap])
  for(i in unique(nchar(hap_comb_array))){
    input_tmp<-which(nchar(hap_comb_array)==i)  
    nmeth<-length(grep("C",hap_comb_array[input_tmp]))
    nunmeth<-length(grep("T",hap_comb_array[input_tmp]))
    mhl_value<-c(mhl_value,i*nmeth/(nunmeth+nmeth))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}

mf<-function(hapinfo){
  mf_value<-c()
  hap_comb_array<-unlist(hapinfo)
  meth_hap=unlist(lapply(hap_comb_array,function(x) unlist(strsplit(x,""))))
  nmeth<-length(grep("C",meth_hap))
  nunmeth<-length(grep("T",meth_hap))
  mf_value<-nmeth/(nunmeth+nmeth)
  return(mf_value)
}

# Figuure 2A
hapinfo<-c(rep("TTTT",16))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2B
hapinfo<-c(rep("CCCC",16))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2C
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
MHL<-mhl(hap_comb(unique(hapinfo)))
MF<-mf(hapinfo)
me<-methentropy(hapinfo)
# Figuure 2D
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2E
hapinfo<-c(rep("CCTT",8),rep("TTCC",8))
MHL<-mhl(hap_comb(unique(hapinfo)))

## Figure 3
gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  GSI<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,mean)))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(data[i,][which(index==group[j])]))/10^(mean(data[i,][which(index==gmax)])))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    GSI<-c(GSI,gsit)
    #  print(c(gmax,gsit))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi)
  return(rlt)
}
ggbarplot<-function(data){
  library("ggplot2")
  newdata<-data.frame(value=as.numeric(data.matrix(data)),Group=rep(colnames(data),each=nrow(data)))
  myData <- aggregate(newdata$value,by =list(type=newdata$Group),
                      FUN = function(x) c(mean = mean(x,na.rm=T),
                                          sd = sd(x,na.rm=T),
                                          sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                          min=min(x,na.rm=T),
                                          max=max(x,na.rm=T),
                                          median=median(x,na.rm=T),
                                          me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
  )
  myData <- do.call(data.frame, myData)
  colnames(myData)=c("type","mean","sd","sem","min","max","median","me")
  ggplot<-ggplot(myData, aes(x =type, y = mean)) +
    geom_bar(position = position_dodge(), stat="identity", fill="blue",width=0.8) +
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.0) +
    theme_bw() +
    theme(panel.grid.major = element_blank())+
    coord_flip()+
    theme(axis.text=element_text(size=10),axis.title=element_text(),axis.text.y = element_text(hjust=0))
  return(ggplot)
}

ggboxplot<-function(data){
  library("ggplot2")
  newdata<-data.frame(value=as.numeric(data.matrix(data)),Group=rep(colnames(data),each=nrow(data)))
  ggplot<-ggplot(newdata, aes(factor(Group),value)) +
    geom_boxplot(aes(fill = factor(Group)),outlier.shape=NA)+
    coord_flip()+
    theme(axis.text=element_text(size=10),axis.title=element_text(),axis.text.y = element_text(hjust=0),legend.position="none")
  return(ggplot)
}

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    dat<-data[,-NaCOL]
  }else{
    dat<-data;
  }
  dat
}   

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

ttestFunction<-function(data,x1,x2){
  data=data+matrix(rnorm(length(data),0.000001,0.000001),nrow(data),ncol(data))
  output<-matrix(NA,dim(data)[1],4)
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),length(na.omit(data[i,x1]))<2,length(na.omit(data[i,x2]))<2,all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]), na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-as.numeric((mean(as.numeric(data[i,x1]),na.rm=T)-mean(as.numeric(data[i,x2]),na.rm=T)))
      output[i,3]<-mean(as.numeric(data[i,x1]),na.rm=T)
      output[i,4]<-mean(as.numeric(data[i,x2]),na.rm=T)
      # print(i)
    }
  }
  
  rownames(output)=rownames(data)
  P.Adj<-p.adjust(output[,1],method="fdr")
  out<-data.frame(output[,1],P.Adj,output[,2:4])
  out<-na.omit(out)
  colnames(out)=c("Pvalue","FDR","Delta","MG1","MG2")
  return(out)
}

f3<-data[mostVarRegion,]
d <- dist(f3, method = "euclidean") # select different distance 
fit <- hclust(d, method="average")  # select different cluster method 
pdf(histFfile)
plot(fit,cex=0.5,hang=-1) # display dendogram
dev.off()
mydata1<-f3
hclustfunc <- function(x) hclust(x, method="ward.D")
distfunc <- function(x) dist(x, method="euclidean")
# perform clustering on rows and columns
cl.row <- hclustfunc(distfunc(mydata1))
cl.col <- hclustfunc(distfunc(t(mydata1)))
gr.row <- cutree(cl.row, 6)
col1 <- brewer.pal(6, "Set1")     # the maximum of brewer.pal is 12
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
col2<-tol21rainbow[as.numeric(as.factor(cl.col$labels))]
library("grDevices")
library("gplots")
col=colorRampPalette(c("blue", "yellow"))(20) 
pdf("heatmap.wgbs.mhl.gsi0.6.without.H1.wbc.cancer.pdf")
heatmap.2(as.matrix(mydata1),hclustfun=hclustfunc, distfun=distfunc,
          RowSideColors=col1[gr.row], 
          ColSideColors=col2,
          Colv=T,Rowv=T,key=T,
          col=col,
          cexCol=0.65,
          keysize=1,
          labRow=NA,
          trace="none",
          density.info="none")
dev.off()


## Figure 4
DataSummary<-function(data){
  rlt<-list()
  NewData<-data.frame(value=as.numeric(data),Group=rep(colnames(data),each=nrow(data)))
  ttest<-t.test(NewData[grep("CP",NewData$Group),1],NewData[grep("Kun",NewData$Group),1],na.rm=T)
  myData <- aggregate(NewData$value,by =list(type=NewData$Group),
                      FUN = function(x) c(mean = mean(x,na.rm=T), 
                                          sd = sd(x,na.rm=T),
                                          sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                          me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
  )
  myData <- do.call(data.frame, myData)
  colnames(myData)=c("type","mean","sd","sem","me")
  rlt$ttest<-ttest
  rlt$myData<-myData
  return(rlt)
}

barplotplasma<-function(myData){
  library("ggplot2")
  # pdf("Lung.cancer.mhl.plasma.groups.107.pdf")
  ylab="Average MHL"
  xlab="Group"
  title="Average MHL in different Groups"
  myData <- within(myData,type <- factor(type,levels=c("NP-Kun","WB","ONT","NLT","LCT","LCP")))
  library("ggplot2")
  ggplot(myData, aes(x =type, y = mean)) +  
    geom_bar(position = position_dodge(), stat="identity", fill="blue",width=0.9) + 
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),size=1.0) +
    ggtitle(title) + 
    theme_bw() +
    theme(panel.grid.major = element_blank())+
    xlab(xlab) +
    ylab(ylab)+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10),axis.text.y = element_text(hjust=0))
}

findConcentration<-function(x,y){
  rlt<-round(rnorm(1,35,2)); # detection limitation (0.001%)
  for(i in 2:length(y)){
    if(x<y[i-1] & x>y[i]){
      rlt<-i
      break
    }
  }
  return(rlt)
}

HeatMap<-function(data){
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  ColSideColors=sidecol(colnames(data))
  heatmap.2(data,trace="none",cexRow = 0.1,cexCol = 0.7, ColSideColors=ColSideColors,density.info="none",col=colors,Colv=F,keysize=0.9, margins = c(5, 10))
}

setwd("C:\\Users\\monod\\MHL")

## Miss Value Statistic
sum(is.na(data))/(nrow(data)*ncol(data))
na<-apply(data,2,function(x) sum(is.na(x)))
sort(na)
barplot(sort(na),horiz=T,las=1,cex.names=0.25,col="blue")
write.table(sort(na),file="missing.value.table.txt",sep="\t")

## re-group/collect the samples #colon cancer
Group<-colnames(data)
# Group[grep("CTR|NC-|Pregn",Group)]<-"NP"
Group[grep("NC-",Group)]<-"NP-Kun"
Group[grep("CTR|Pregn",Group)]<-"NP-Dennis"
Group[grep("6-P-",Group)]<-"CCP"
Group[grep("6-T-|HCT116|Colon_Tumor_Primary|CTT-|metastasis_colon|tumor_colon",Group)]<-"CCT"
Group[grep("normal_colon|N37-Colon|SG-01",Group)]<-"NCT"
# Lung cancer
Group[grep("7-P-",Group)]<-"LCP"
Group[grep("7-T-|adenocarcinoma_lung|tumor_lung",Group)]<-"LCT"
Group[grep("LG-01|N37-Lung|normal_lung",Group)]<-"NLT"
# Pancreatic cancer
Group[grep("PC-P-",Group)]<-"PCP"
Group[grep("PC-T-",Group)]<-"PCT"
Group[grep("PA-01|N37-Pancreas",Group)]<-"NPT"
# WB and Normal tissues
Group[grep("WB-",Group)]<-"WB"
Group[grep("STL|N37-|methylC-",Group)]<-"ONT"

# data proxy
newdata<-data
colnames(newdata)<-Group
newdata<-data.frame(newdata,check.names=F)

## action
YY<-grep("NP-Kun|CCP|CCT|NCT|LCP|LCT|NLT|PCP|PCT|NPT|WB|ONT",Group)
Newdata<-newdata[,YY]
Data_tmp1<-data[,YY]   # heatmap plot with ID

# validate Newdata is same as Data_tmp1
# Newdata[1:10,1:10]
# Data_tmp1[1:10,1:10]
# dim(Newdata)
# dim(Data_tmp1)

idx<-sapply(colnames(Newdata),function(x) unlist(strsplit(x,"[.]"))[1])
colnames(Newdata)<-idx
Newdata<-data.matrix(Newdata)
head(Newdata)
colnames(Newdata)
table(colnames(Newdata))
dim(Newdata)

x1<-grep("NP-Kun",colnames(Newdata))
x2<-grep("CCP",colnames(Newdata))
x3<-grep("CCT",colnames(Newdata))
x4<-grep("NCT",colnames(Newdata))
x5<-grep("LCP",colnames(Newdata))
x6<-grep("LCT",colnames(Newdata))
x7<-grep("NLT",colnames(Newdata))
x8<-grep("PCP",colnames(Newdata))
x9<-grep("PCT",colnames(Newdata))
x92<-grep("NPT",colnames(Newdata))
x10<-grep("WB",colnames(Newdata))
x11<-grep("ONT",colnames(Newdata))

# Data_tmp1<-Data_tmp1[,c(x1,x2,x3,x4,x10,x11)]    # heatmap plot with ID
# Newdata<-Newdata[,c(x1,x2,x3,x4,x10,x11)]        # for colon cancer
Newdata<-Newdata[,c(x1,x5,x6,x7,x10,x11)]          # for lung cancer
# Newdata<-Newdata[,c(x1,x8,x9,x92,x10,x11)]         # for pancreatic cancer

table(colnames(Newdata))
# remove outliner of normal tissues
# Tmp<-Newdata[,grep("NCT|ONT",colnames(Newdata))]
# Tmp[Tmp>0.5]<-NA
# Newdata[,grep("NCT|ONT",colnames(Newdata))]<-Tmp

# validate Newdata is same as Data_tmp1
# Newdata[1:10,1:10]
# Data_tmp1[1:10,1:10]
# dim(Newdata)
# dim(Data_tmp1)
colnames(Newdata)
table(colnames(Newdata))
dim(Newdata)
colnames(Data_tmp1)

thres1<-0.05
thres2<-0.5
target1.colon<-which(apply(Newdata,1,function(x)  mean(x[grep("WB",colnames(Newdata))],na.rm=T)< thres1 && mean(x[grep("LCT",colnames(Newdata))],na.rm=T) > thres2 )) 
# Cancer tissue > 0 and normal plasma =0
length(target1.colon)
Newdata2=Newdata[target1.colon,]

# head(Newdata2)
# head(Newdata2_tmp1)
heatmap<-HeatMap(Newdata2)

rlt<-DataSummary(Newdata2)
myData<-rlt$myData
# myData <- myData[-which(myData$type=="CCT"),]
myData
barplotplasma(myData)

# hclust and the sub-genomic region biomarkers(row is gene region)
dist<-dist((Newdata2))
hclust<-hclust(dist)
plot(hclust,cex=0.5)

k=80
cut<-cutree(hclust,k=k)
xx<-cut[match(rownames(Newdata2)[heatmap$rowInd],names(cut))]
yy<-unique(xx)
# sample size in each subgroup
zz<-table(xx)[match(yy,names(table(xx)))]

# Group 2: 
Newdata2<-Newdata2[match(names(xx[which(xx %in% c(15))]),rownames(Newdata2)),]
input.long<-data.frame(variable=names(colMeans(Newdata2,na.rm = T)),value=colMeans(Newdata2,na.rm = T))
table(input.long$variable)
input.long <- within(input.long,variable <- factor(variable,levels=c("NP-Kun","WB","ONT","NLT","LCT","LCP")))
ggplot(aes(y = value, x = variable), data = input.long) + geom_boxplot(outlier.shape =16,outlier.colour="red")+theme_bw() + ylim(0,1)+ theme(
  plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank()
)

## Figure 5
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

#################################################
################## CCP ##########################
#################################################
par(mfrow=c(2,2))
hist(CCP[2,],breaks=30,xlim=c(0,60),col="darkblue",border="darkblue",xlab="Colon",ylab="Counts of cs-MHL",main="Colon")
hist(CCP[11,],breaks=10,xlim=c(0,60),col="darkblue",border="darkblue",xlab="CT",ylab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)

Zmax<-matrix(nrow=nrow(CCP),ncol=ncol(CCP))
for(i in 1:ncol(CCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,60)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    zmp <- (CCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(CCP)
colnames(Zmax)=colnames(CCP)
Zmax

par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i])
}

par(mfrow=c(2,2))
library("beeswarm")
input<-list(CCP=Zmax[2,],WB=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:4,pch=16,method="center",ylab="Z-Score")


p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="bonferroni"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)

sum(p1[2,]>-log(0.05,10))
length(unique(c(which(p1[2,]>-log(0.05,10)),which(p1[11,]>-log(0.05,10)))))

library("gplots")
heatmap.2(p1,Rowv=F,Colv=F,density.info="none",trace="none")


write.table(p1,file="colon.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)
getwd()

# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==2))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))/ncol(Zmax)

### Figure 
par(mfrow=c(4,3),mar=c(4,4,2,1))
AUC<-c()
for(k in 1:11){
  xx<-c()
  sen<-c()
  spe<-c()
  for(fdr in seq(0,1,by=0.01)){
    xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
    sen<-c(sen,xx/30)
    spe<-c(spe,1-fdr)
  }
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main=rownames(Zmax)[k])
  AUC<-c(AUC,mean(sen))
}

Zmax[k,]<-Zmax[2,]+Zmax[11,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Colon+CT")

# barplot of AUC for CRC cancer
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Colon+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")


#################################################
################## LCP ##########################
#################################################

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\tissue-of-origin-mapping")
normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))
Zmax<-matrix(nrow=nrow(LCP),ncol=ncol(LCP))
for(i in 1:ncol(LCP)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,29)
    idx2<-idx[which(! idx %in% idx1)]
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    zmp <- (LCP[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}
rownames(Zmax)=rownames(LCP)
colnames(Zmax)=colnames(LCP)
Zmax
par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i])
}

library("beeswarm")
input<-list(LCP=Zmax[6,],WB=Zmax[10,],CT=Zmax[11,])
beeswarm(input,col = 2:4,pch=16,method="center",ylab="Z-Score")
par(mfrow=c(2,2))
hist(LCP[6,],breaks=20,xlim=c(0,60),col="darkblue",border="darkblue",ylab="Counts of samples",xlab="Counts of ls-MHL",main="Colon")
hist(LCP[11,],breaks=8,xlim=c(0,60),col="darkblue",border="darkblue",ylab="Counts of samples",xlab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)
# bar plot show Z-score for one lung sample
barplot(Zmax[,2],col="darkblue",ylab="Z-score")
file=list.file("*txt")
for(i in file){
  d<-read.table(i)
}
# merge difference cancers into one cancer group and take all the notehrs.
sort(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==2)))
length(which(apply(Zmax,2,which.max)==6))/ncol(Zmax)
length(which(apply(Zmax,2,which.max)==11))/ncol(Zmax)
length(c(which(apply(Zmax,2,which.max)==11),which(apply(Zmax,2,which.max)==6)))/ncol(Zmax)
p<-pnorm(-abs(Zmax))
p1<-matrix(p.adjust(p,method="bonferroni"),nrow=11,byrow=F)
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1<--log(p1,10)
length(unique(c(which(p1[6,]>-log(0.001,10)),which(p1[11,]>-log(0.001,10)))))
write.table(p1,file="lung.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)
## Figure C
par(mfrow=c(4,3),mar=c(4,4,2,1))
AUC<-c()
for(k in 1:11){
  xx<-c()
  sen<-c()
  spe<-c()
  for(fdr in seq(0,1,by=0.01)){
    xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
    sen<-c(sen,xx/30)
    spe<-c(spe,1-fdr)
  }
  AUC<-c(AUC,mean(sen))
  plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main=rownames(Zmax)[k])
}
Zmax[k,]<-Zmax[6,]+Zmax[11,]
xx<-c()
sen<-c()
spe<-c()
for(fdr in seq(0,1,by=0.01)){
  xx<-sum(Zmax[k,]> quantile((normal[k,]-mean(normal[k,]))/(sd(normal[k,])/sqrt(length(normal[k,]))),fdr))
  sen<-c(sen,xx/30)
  spe<-c(spe,1-fdr)
}
plot(sen~1-spe,col="red",pch=16,cex=1,xlab="1-specificity",ylab="sensitivity",main="Lung+CT")

# barplot of AUC for Lung cancer
par(mfrow=c(1,1))
AUC<-c(AUC,mean(sen))
names(AUC)<-c(rownames(Zmax),"Lung+CT")
xx<-barplot(AUC,col="blue",ylim=c(0,1))
text(x = xx, y = AUC, label = round(AUC,3), pos = 3, cex = 0.8, col = "red")


####################################################
################## Normal ##########################
####################################################

normal<-data.matrix(read.table("normal-ref.txt",sep="\t",head=T,row.names=1,as.is=T))
CCP<-data.matrix(read.table("CCP.txt",sep="\t",head=T,row.names=1,as.is=T))
LCP<-data.matrix(read.table("LCP.txt",sep="\t",head=T,row.names=1,as.is=T))

par(mfrow=c(3,4))
for(i in 1:(nrow(normal))){
  hist(normal[i,],col="darkblue",breaks=30,xlim=c(0,50),xlab=rownames(normal)[i],border ="darkblue",main="",ylab="Counts of rsMHL",cex.axis=1.15,cex.lab=1.5)  
  abline(v=30,lty=3,col="red",lwd=2)
}


Zmax<-matrix(nrow=nrow(normal),ncol=ncol(normal))
for(i in 1:ncol(normal)){
  for(k in 1:nrow(normal)){
    idx<-1:ncol(normal)
    idx1<-sample(idx,75)
    idx2<-idx[which(! idx %in% idx1)]
    
    Mean<-mean(normal[k, idx1])
    SD<-sd(normal[k, idx1])
    
    Z<-c()
    for(p in normal[k,idx2]){
      z <- (p - Mean)/(SD/sqrt(length(idx1)))
      Z<-c(Z,z)
    }
    
    zmp <- (normal[k,i] - mean(normal[k, idx]))/(sd(normal[k, idx])*sqrt((length(normal[k, idx])-1)/(length(normal[k, idx]))))
    Zmax[k,i]=zmp 
  }
}

rownames(Zmax)=rownames(normal)
colnames(Zmax)=colnames(normal)
Zmax



p<-pnorm(-abs(Zmax))
rownames(p1)<-rownames(p)
colnames(p1)<-colnames(p)
p1
p1<--log(p1,10)

sum(p[2,]<0.05)
sum(p[6,]<0.05)

write.table(p1,file="normal.pvalue.matrix.txt",col.names=NA,row.names=T,sep="\t",quote=F)


##############################################################
################## CCP.Z-score.txt ##########################
##############################################################

Zmax<-read.table("CCP.Z-score.txt")
par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i],col = "darkblue")
}

Zmax<-read.table("LCP.Z-score.txt")
par(mfrow=c(5,6),mar=c(2,2,2,2))
for(i in 1:ncol(Zmax)){
  barplot(Zmax[,i],main=colnames(Zmax)[i],col = "darkblue")
}

