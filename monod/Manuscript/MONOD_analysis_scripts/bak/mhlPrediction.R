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

Zscore<-function(ccpr,npr){
  Zmax<-matrix(0,nrow=nrow(ccpr),ncol=ncol(ccpr))
  for(i in 1:ncol(ccpr)){
    for(k in 1:nrow(npr)){
      idx<-1:ncol(npr)
      zmp <- (ccpr[k,i] - mean(npr[k, idx]))/(sd(npr[k, idx])*sqrt((length(npr[k, idx])-1)/(length(npr[k, idx]))))
      Zmax[k,i]=zmp 
    }
  }
  rownames(Zmax)=rownames(ccpr)
  colnames(Zmax)=colnames(ccpr)
  Zmax
}

rename2<-function(data){
  saminfo<-read.table("/media/NAS1/shg047/monod/hapinfo/N37Salk.saminfo",sep="\t")
  head(saminfo)
  colnames(data)<-gsub(".allchrs.rmdup.bam.mhbs","",colnames(data))
  colnames(data)<-gsub(".allchrs.sorted.clipped.bam.mhbs","",colnames(data))
  colnames(data)[grep("age|new|centenarian",colnames(data))]<-"WBC"
  colnames(data)[grep("X7.T",colnames(data))]<-"LCT"
  colnames(data)[grep("X6.T",colnames(data))]<-"CCT"
  colnames(data)[grep("methylC-seq_h1",colnames(data))]<-"H1"
  colnames(data)[grep("normal_lung",colnames(data))]<-"Lung"
  colnames(data)[grep("normal_prostate",colnames(data))]<-"Prostate"
  colnames(data)[grep("normal_brain",colnames(data))]<-"Brain"
  colnames(data)[grep("normal_colon|Colon_Primary_Normal",colnames(data))]<-"Colon"
  colnames(data)[grep("normal_breast",colnames(data))]<-"Breast"
  colnames(data)[grep("normal_liver",colnames(data))]<-"Liver"
  colnames(data)[grep("normal_CD19",colnames(data))]<-"WBC"
  colnames(data)[grep("Frontal_cortex_normal",colnames(data))]<-"Brain"
  colnames(data)[grep("fetal_heart_cells",colnames(data))]<-"Heart"
  colnames(data)[grep("SRX381710_normal_placenta|fPlacenta_cells",colnames(data))]<-"Placenta"
  colnames(data)[grep("adult_CD*",colnames(data))]<-"WBC"
  colnames(data)[grep("fAdrenal_cells",colnames(data))]<-"Adrenal"
  colnames(data)[grep("fThymus_cells",colnames(data))]<-"Thymus"
  colnames(data)[grep("fMuscle_Leg__cells|Muscle_Trunk__cells",colnames(data))]<-"Muscle"
  colnames(data)[grep("fThymus_cells",colnames(data))]<-"Thymus"
  colnames(data)[grep("fStomach_cells",colnames(data))]<-"Stomach"
  colnames(data)[grep("N37",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("N37",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("STL",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("STL",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("Indx",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("Indx",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("X6.P|RRBS.6P",colnames(data))]<-"CCP"
  colnames(data)[grep("X7.P|RRBS.7P",colnames(data))]<-"LCP"
  colnames(data)[grep("X7.P|RRBS.7P",colnames(data))]<-"LCP"
  colnames(data)[grep("BS-seq-P1-N|BS-seq-P2-N",colnames(data))]<-"Kidney"
  colnames(data)[grep("BS-seq-P1-T|BS-seq-P2-T",colnames(data))]<-"PKIRCT"
  colnames(data)[grep("tumor_liver",colnames(data))]<-"PLIHCT"
  colnames(data)[grep("tumor_glioblastoma|tumor_neuroblastoma",colnames(data))]<-"PGBMT"
  colnames(data)[grep("CD19_cells",colnames(data))]<-"WBC"
  colnames(data)[grep("Colon_Tumor_Primary|SRX381569_tumor_colon",colnames(data))]<-"PCOADT"
  colnames(data)[grep("SRX381621_tumor_breast",colnames(data))]<-"PBRCAT"
  colnames(data)[grep("tumor_prostate",colnames(data))]<-"PPRADT"
  colnames(data)[grep("SRX381722_small_cell_tumor_lung|SRX381719_squamous_cell_tumor_lung|SRX381716_adenocarcinoma_lung",colnames(data))]<-"PLCT"
  colnames(data)[grep("muscle",colnames(data))]<-"Muscle"
  colnames(data)[grep("Frontal_cortex_AD",colnames(data))]<-"FrontalAD"
  colnames(data)[grep("intestine",colnames(data))]<-"Intestine"
  colnames(data)[grep("H1",colnames(data))]<-"H1"
  return(data)
}


prediction<-function(data,bio,tt=0.6){
  rlt<-list()
  input<-data[na.omit(match(rownames(bio),rownames(data))),]
  bio<-bio[na.omit(match(rownames(input),rownames(bio))),]
  factors<-c("Brain","Colon","Intestine","Kidney","Liver","Lung","Pancreas","Spleen","Stomach","WBC")
  testcounts<-apply(input,2,function(x) table(factor(unlist(lapply(bio[match(rownames(input)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,",")))),levels=factors)))
  backcounts<-table(unlist(lapply(bio[match(rownames(input),rownames(bio)),]$group,function(x) unlist(strsplit(x,",")))))
  prediction<-apply(apply(testcounts,2,function(x) x/backcounts),2,function(x) rownames(backcounts)[which.max(x)])
  score<-apply(testcounts,2,function(x) x/backcounts)
  rlt$testcounts=testcounts
  rlt$backcounts=backcounts
  rlt$prediction=prediction
  rlt$score=score
  return(rlt)
}


gsi<-function(matrix){
  narow<-which(unlist(apply(matrix,1,function(x) all(is.na(x)))))
  if(length(narow)>0){
    matrix<-matrix[-narow,]
  }
  group = names(table(colnames(matrix)))
  index = colnames(matrix)
  GSI <- c()
  gmaxgroup <- c()
  refMean<-c()
  refMax<-c()
  for (i in 1:nrow(matrix)) {
    gsit <- 0    
    refmean<-tapply(as.numeric(matrix[i, ]), index, function(x) mean(x, na.rm = T))
    refmax<-refmean[which.max(refmean)]
    gmax <- names(refmax)
    for (j in 1:length(group)) {
      tmp <- (1 - 10^(mean(na.omit(as.numeric(matrix[i, which(index == group[j])])), na.rm = T))/10^(mean(na.omit(as.numeric(matrix[i,which(index == gmax)])))))/(length(group) - 1)
      gsit <- c(gsit, tmp)
    }
    if(sum(refmean>0.3,na.rm=T)>1){
      gmax<-gsub(" ","",paste(unique(c(gmax,names(refmean)[which(refmean>0.3)])),',',collapse =""))
    }
    gmaxgroup <- c(gmaxgroup, gmax)
    GSI <- c(GSI, sum(gsit, na.rm = T))
    refMean<-c(refMean,gsub(" ","",paste(round(refmean,5),',',collapse ="")))
    refMax<-c(refMax,refmax)    
  }  
  rlt = data.frame(region = rownames(matrix), group = gmaxgroup, GSI = GSI, refMax=refMax,refMean=refMean)
  return(rlt)
}



HeatMap<-function(data,cexRow=0.5,cexCol=0.5){
  #  note: this function include correlation based heatmap (pearson or spearman)
  #  data: row is gene and column is sample
  #  colname and rowname cannot be NULL  
  #  Usage example:
  #  test<- matrix(runif(100),nrow=20)
  #  colnames(test)=c("A","A","A","B","B")
  #  rownames(test)=paste("Gene",1:20,sep="")
  #  HeatMap(test)
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  
  Hclust=function(x){hclust(x,method="complete")}
  Distfun=function(x){as.dist((1 - cor(t(x),method = "spearman")))}
  
  ColSideColors=sidecol(colnames(data))
  
  heatmap.2(data,trace="none", 
            hclust=Hclust,
            distfun=Distfun, 
            cexRow = cexRow, cexCol = cexCol,
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=T,Rowv = TRUE,
            keysize=0.9, margins = c(5, 10)
  )
}

ZscorePredictionPlasma<-function (test, ncp, bio, tt = 0.001)
{
  rlt <- list()
  test <- test[rownames(test) %in% rownames(bio), ]
  ncp <- ncp[rownames(ncp) %in% rownames(bio), ]
  ccp <- test[rownames(test) %in% rownames(ncp), ]
  
  # test1, fragment vs background
  bio<-bio[na.omit(match(rownames(test),rownames(bio))),]
  backcounts<-table(unlist(lapply(bio[match(rownames(test),rownames(bio)),]$group,function(x) unlist(strsplit(x,",")))))
  factors<-c("Brain","Colon","Intestine","Kidney","Liver","Lung","Pancreas","Stomach")
  testcounts<-apply(test,2,function(x) table(factor(unlist(lapply(bio[match(rownames(test)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,",")))),levels=factors)))
  score<-apply(testcounts,2,function(x) x/backcounts)
  prediction<-apply(score,2,function(x) rownames(backcounts)[which.max(x)])
  
  # test1, cancer plasa  vs normal plasma
  npr <- apply(ncp, 2, function(x) table(unlist(lapply(bio[match(rownames(ncp)[which(x >
                                                                                       tt)], rownames(bio)), ]$group, function(x) unlist(strsplit(x,
                                                                                                                                                  ","))))))
  ccpr <- apply(ccp, 2, function(x) table(unlist(lapply(bio[match(rownames(ccp)[which(x >
                                                                                        tt)], rownames(bio)), ]$group, function(x) unlist(strsplit(x,
                                                                                                                                                   ","))))))
  # output
  rlt$testcounts<-testcounts
  rlt$backcounts<-backcounts
  rlt$score<-score
  rlt$npn<-npr
  rlt$ccpn<-ccpr
  rlt$nprscore <- Zscore(npr, npr)
  rlt$ccprscore <- Zscore(ccpr, npr)
  rlt1 <- apply(Zscore(ccpr, npr), 2, function(x) rownames(npr)[which.max(x)])
  rlt2 <- apply(Zscore(npr, npr), 2, function(x) rownames(npr)[which.max(x)])
  rlt$test <- rlt1
  rlt$background <- rlt2
  return(rlt)
}



data<-read.table(opt$file,head=T,row.names=T,sep="\t")

## MHL
#data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.mhl.mhbs",head=T,row.names=1,sep="\t",check.names = F)
data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.mhl.mhbs.Guo.txt",head=T,row.names=1,sep="\t",check.names = F)
#rownames(data)<-unlist(lapply(rownames(data),function(x) paste(unlist(strsplit(x,"[:|-]"))[1],":",unlist(strsplit(x,"[:|-]"))[2],"-",unlist(strsplit(x,"[:|-]"))[3],sep="")))
#write.table(data,file="/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.mhl.mhbs.Guo.txt",col.names=NA,row.names=T,sep="\t",quote=F)

## AMF
#data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.amf.mhbs",head=T,row.names=1,sep="\t",check.names = F)
#data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.amf.mhbs.Guo.txt",head=T,row.names=1,sep="\t",check.names = F)
#rownames(data)<-unlist(lapply(rownames(data),function(x) paste(unlist(strsplit(x,"[:|-]"))[1],":",unlist(strsplit(x,"[:|-]"))[2],"-",unlist(strsplit(x,"[:|-]"))[3],sep="")))
#write.table(data,file="/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.amf.mhbs.Guo.txt",col.names=NA,row.names=T,sep="\t",quote=F)


# rename and select tissues to run gsi and then get tissue specific biomarkers
data<-rename2(data)
table(colnames(data))
#Data<-data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|Heart|WBC|CCP|LCP|NP",colnames(data))] 
Data<-data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(data))] 
colnames(Data)
colnames(Data)<-unlist(lapply(colnames(Data),function(x) unlist(strsplit(x,"[.]"))[1]))
#gsirlt<-gsi(Data)
#save(gsirlt,file="gsirlt.mhl.may30.RData")

# subselect better tissue specifci makres. We find with the following parameters, tissue have highest prediction accuracy.
# When I found there are huge number of markers shared by WBC, lung and Colon, I am thinking why not remove WBC makres and then check non-WBC signals in plasma
# Another thing wwe need to do is to remove one Lung and colon samples which are quite similiar with WBC in the cluster analysis, maybe they are labelled mistakely.  
load("gsirlt.mhl.may30.RData")
gsirlt<-gsirlt[-which(unlist(lapply(gsirlt$group,function(x) length(grep("WBC|Spleen",x))>0))),]
gsirlt1<-subset(gsirlt,refMax>0.3 & GSI>0.6)  # refMax=c(0.3-0.6), GSI=c(0.5-0.7)
gsirlt1$group<-as.character(gsirlt1$group)
bio<-gsi2bio(gsirlt1)
rlt1<-prediction(Data,bio,tt=0.3)              # tissue=c(0.3-0.6)
rlt2<-prediction(data,bio,tt=0.3)              # tissue=c(0.3-0.6)
rlt3<-prediction(rename2(data),bio,tt=0.3)      # reanmed data
sum(rlt1$prediction==names(rlt1$prediction))/length(rlt1$prediction)

# colon, lung and WBC are shared large number common hyermethyalted regions compared with othter tissues.
pdf("heatmap-tovb.pdf")                      # heatmap to all the samples with occurred tissue marker vs background
heatmaprlt<-HeatMap(rlt2$score)
heatmaprlt<-HeatMap(t(rlt2$score))
dev.off()

# prepare RRBS dataset
#rrbsdata<-read.table("/media/NAS1/shg047/monod/mhl/0530/RRBS.getHaplo.mhl.mhbs",head=T,row.names=1,sep="\t",check.names = F)
#rownames(rrbsdata)<-unlist(lapply(rownames(rrbsdata),function(x) paste(unlist(strsplit(x,"[:|-]"))[1],":",unlist(strsplit(x,"[:|-]"))[2],"-",unlist(strsplit(x,"[:|-]"))[3],sep="")))
#write.table(rrbsdata,file="/media/NAS1/shg047/monod/mhl/0530/RRBS.getHaplo.mhl.mhbs.Guo.txt",col.names=NA,row.names=T,sep="\t",quote=F)
rrbsdata<-read.table("/media/NAS1/shg047/monod/mhl/0530/RRBS.getHaplo.mhl.mhbs.Guo.txt",head=T,row.names=1,sep="\t",check.names = F)

# phase I predict plasma with model of solid tissue, then all the plasma samples should be category to WBC
rrbsdata<-rrbsdata[,grep("CRC-P|LC-P|NC-P",colnames(rrbsdata))]
rlt1<-prediction(rrbsdata,bio,tt=0.6)             

# heatmap to show the tissue-specific markers occured in plasmas(all plasma)
pdf("plasma.tsMHL-counts-in-plasma.pdf")
HeatMap(rlt1$testcounts)
HeatMap(t(rlt1$testcounts),cexRow=0.25,cexCol=0.5)
HeatMap(rlt1$score)
HeatMap(t(rlt1$score),cexRow=0.25,cexCol=0.5)
dev.off()

# separate data into train and test for normal plasma
np<-rrbsdata[,grep("NC-P|NCP|NP",colnames(rrbsdata))]  
train<-sample(1:75,45)
ncp<-np[,train]
nncp<-np[,(1:75)[-train]]

ccp<-rrbsdata[,unique(grep(".6P|X6.P|CCP|CRC-P",colnames(rrbsdata)))]  
lcp<-rrbsdata[,unique(grep(".7P|X7.P|LC-P",colnames(rrbsdata)))] 

for(i in seq(0,1,by=0.001)){
rlt1<-ZscorePredictionPlasma(nncp,ncp,bio,tt=i)   # tt is threshold to binary the MHL matrix, lcp: colon cancer plasma, ncp: normal plasma trainning dataset
rlt2<-ZscorePredictionPlasma(ccp,ncp,bio,tt=i)   # tt is threshold to binary the MHL matrix, ccp: colon cancer plasma, ncp: normal plasma trainning dataset
rlt3<-ZscorePredictionPlasma(lcp,ncp,bio,tt=i)   # tt is threshold to binary the MHL matrix, lcp: colon cancer plasma, ncp: normal plasma trainning dataset

prin1<-sum(unlist(apply(rlt1$score,2,function(x) rownames(rlt1$score)[which.max(x)]))=="Colon")
prin2<-sum(unlist(apply(rlt2$score,2,function(x) rownames(rlt1$score)[which.max(x)]))=="Lung")
print(c(i,prin1,prin2))
}



