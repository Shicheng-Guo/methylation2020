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
  testcounts<-apply(input,2,function(x) table(unlist(lapply(bio[match(rownames(input)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,","))))))
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

HeatMap<-function(data){
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
            cexRow = 0.5, cexCol = 0.5,
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=T,Rowv = TRUE,
            keysize=0.9, margins = c(5, 10)
  )
}

TSNEAnalysis<-function(mydata,phen,prefix="TSNE"){
  print("Missing Value: NA is not allowed in the input matrix!")
  print("Please be sure row is individual, column is variable!")
  print("t-sne analysis: N rows (objects) x P columns (variables) [same as PCA], be sure rownames should be unique")
  library("tsne")
  output=paste(prefix,".tsne-1.pdf",sep="")
  pdf(output)
  colors = rainbow(length(unique(phen)))
  names(colors) = unique(phen)
  ecb = function(x,y){ plot(x,t='n',xlab="Coordinate 1",ylab="Coordinate 2"); text(x,labels=phen, col=colors[phen]) }
  tsne_iris = data.frame(tsne(mydata, epoch_callback = ecb, perplexity=50))
  dev.off()
  
  colnames(tsne_iris)<-c("xtsne","ytsne")
  output=paste(prefix,".tsne.pdf",sep="")
  pdf(output)
  chart = ggplot(data.frame(tsne_rlt), aes(xtsne,ytsne))+geom_point(size=1,alpha=1,aes(colour = factor(phen1)))+ggtitle("tSNE dimensions colored by digit")
  # change the size and alpha could make great beautful figure
  chart
  dev.off()
  
}

PCAPlot<-function(data,pheno,output,multifigure=T){
  print("Please be sure row is individual, column is variable!")
  print("Missing Value: NA is not allowed in the input matrix!")
  nacol<-which(unlist(apply(data,2,function(x) all(is.na(x)))))
  if(length(nacol)>0){
    matrix<-matrix[,-narow]
  }
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
  legend("topright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",cex=0.6)
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("bottomright",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",cex=0.6)
  dev.off()
}

RowNARemove<-function(data,missratio=0.9){
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

ColNARemove<-function(data,missratio=0.9){
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


optGsi<-function(Data,gsirlt){
  Acc<-c()
  for(i in seq(0.1,0.9,by=0.1)){
    for(i in seq(0.1,0.9,by=0.1)){
      for(k in seq(0.1,0.9,by=0.1)){
        gsirlt1<-subset(gsirlt,refMax>0.3 & GSI>0.6)  # refMax=c(0.3-0.6), GSI=c(0.5-0.7)
        gsirlt1$group<-as.character(gsirlt1$group)
        bio<-gsi2bio(gsirlt1)
        rlt1<-prediction(Data,bio,tt=0.3)              # tissue=c(0.3-0.6)
        Acc<-rbind(Acc,c(i,j,k,sum(rlt1$prediction==names(rlt1$prediction))/length(rlt1$prediction)))
      }
    }
  }
}

data<-read.table(opt$file,head=T,row.names=T,sep="\t")

#data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.mhl.mhbs",head=T,row.names=1,sep="\t",check.names = F)
#data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.mhl.mhbs.Guo.txt",head=T,row.names=1,sep="\t",check.names = F)
#rownames(data)<-unlist(lapply(rownames(data),function(x) paste(unlist(strsplit(x,"[:|-]"))[1],":",unlist(strsplit(x,"[:|-]"))[2],"-",unlist(strsplit(x,"[:|-]"))[3],sep="")))
#write.table(data,file="/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.mhl.mhbs.Guo.txt",col.names=NA,row.names=T,sep="\t",quote=F)

data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.amf.mhbs",head=T,row.names=1,sep="\t",check.names = F)
data<-read.table("/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.amf.mhbs.Guo.txt",head=T,row.names=1,sep="\t",check.names = F)
rownames(data)<-unlist(lapply(rownames(data),function(x) paste(unlist(strsplit(x,"[:|-]"))[1],":",unlist(strsplit(x,"[:|-]"))[2],"-",unlist(strsplit(x,"[:|-]"))[3],sep="")))
write.table(data,file="/media/NAS1/shg047/monod/mhl/0530/WGBS.getHaplo.amf.mhbs.Guo.txt",col.names=NA,row.names=T,sep="\t",quote=F)

# rename and select tissues to run gsi and then get tissue specific biomarkers
data<-ColNARemove(data,0.9)
data<-RowNARemove(data,0.9)

data<-rename2(data)
table(colnames(data))
Data<-data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|Heart|WBC|CCP|LCP|NP",colnames(data))] 
colnames(Data)
colnames(Data)<-unlist(lapply(colnames(Data),function(x) unlist(strsplit(x,"[.]"))[1]))

# non-supervisored cluster(anyone think about cluster analysis to binrary methylation dataset?)
# top5000 sd cluster(anyone think about cluster analysis to binrary methylation dataset?)
sd<-unlist(apply(data,1,function(x) sd(na.omit(x))))
topsd5000<-data[order(sd,decreasing=T)[1:5000],]

pdf("heatmap.topsd5000.pdf")
HeatMap(topsd5000)
HeatMap(t(topsd5000))
dev.off()

# top5000 GSI cluster(anyone think about cluster analysis to binrary methylation dataset?)
TSNEAnalysis(data.matrix(topsd5000),phen=colnames(data),prefix="TSNE")
Rscript ~/work/Alice/mouse/scMeth/hapinfo/heatmap.R
# PCA
PCAPlot(t(data),phen=rownames(rename2(data)),output="pca-new-amf-may30.pdf")
gsirlt<-gsi(Data)
save(gsirlt,file="gsirlt.amf.may30.RData")

# subselect better tissue specifci makres. We find with the following parameters, tissue have highest prediction accuracy.
Acctrain<-optGsi(Data,gsirlt)                 # find best refMax,GSI and tt threshold
gsirlt1<-subset(gsirlt,refMax>0.3 & GSI>0.6)  # refMax=c(0.3-0.6), GSI=c(0.5-0.7)
gsirlt1$group<-as.character(gsirlt1$group)
bio<-gsi2bio(gsirlt1)
rlt1<-prediction(Data,bio,tt=0.3)              # tissue=c(0.3-0.6)
sum(rlt1$prediction==names(rlt1$prediction))/length(rlt1$prediction)

rlt2<-prediction(data,bio,tt=0.3)              # tissue=c(0.3-0.6)
sum(rlt2$prediction==names(rlt2$prediction))/length(rlt2$prediction)
rlt3<-prediction(rename2(data),bio,tt=0.3)      # reanmed data
sum(rlt1$prediction==names(rlt1$prediction))/length(rlt1$prediction)
rlt1$predict

# colon, lung and WBC are shared large number common hyermethyalted regions compared with othter tissues.
pdf("heatmap-tovb-amf.pdf")                      # heatmap to all the samples with occurred tissue marker vs background
heatmaprlt<-HeatMap(rlt2$score)
heatmaprlt<-HeatMap(t(rlt2$score))
dev.off()
save.image(file="save.image.RData")

