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

heatmap<-function(mydata,output,cexRow=0.5,cexCol=0.5){
  hclustfunc <- function(x) hclust(x, method="complete")
  distfunc <- function(x) dist(x, method="euclidean")
  # perform clustering on rows and columns
  cl.row <- hclustfunc(distfunc(mydata))
  cl.col <- hclustfunc(distfunc(t(mydata)))
  # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
  gr.row <- cutree(cl.row, 6)
  # gr.col <- cutree(cl.col, 5)
  # require(RColorBrewer)
  col1 <- brewer.pal(6, "Set1")     # the maximum of brewer.pal is 12
  #col2 <- brewer.pal(5, "Pastel1")
  #tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  if(length(unique(colnames(mydata)))>25){
    tol21rainbow<-rainbow(length(unique(colnames(mydata))))
  }else{
    tol21rainbow <- c("dodgerblue2","#E31A1C", # red
                      "green4",
                      "#6A3D9A", # purple
                      "#FF7F00", # orange
                      "black","gold1",
                      "skyblue2","#FB9A99", # lt pink
                      "palegreen2",
                      "#CAB2D6", # lt purple
                      "#FDBF6F", # lt orange
                      "gray70", "khaki2",
                      "maroon","orchid1","deeppink1","blue1","steelblue4",
                      "darkturquoise","green1","yellow4","yellow3",
                      "darkorange4","brown")
  }
  col2<-tol21rainbow[as.numeric(as.factor(cl.col$labels))]
  col=colorRampPalette(c("blue", "yellow"))(20) 
  pdf(output)
  par(mar=c(5,5,5,0))
  heatmaprlt<-heatmap.2(as.matrix(mydata),hclustfun=hclustfunc, distfun=distfunc,
                        RowSideColors=col1[gr.row], 
                        ColSideColors=col2,
                        cexRow=cexRow,cexCol=cexCol,
                        labRow=F,
                        trace="none",
                        col=col,
                        na.color="black",
                        density.info="none")
  legend=unique(data.frame(col2,cl.col$labels))
  legend(x=0.85,y=0.8,legend=legend[,2],col=as.character(legend[,1]),pch=15,cex = 0.5,bty="n")
  dev.off()
  return(heatmaprlt)
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
  library("RColorBrewer")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  
  Hclust=function(x){hclust(x,method="complete")}
  Distfun=function(x){as.dist((1 - cor(t(x),use="pairwise.complete.obs",method = "spearman")))}
  
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

ZscorePredictionPlasma<-function (test, ncp, bio, tt = 0.001){
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
  # test2, cancer plasa  vs normal plasma
  npr <- apply(ncp, 2, function(x) table(unlist(lapply(bio[match(rownames(ncp)[which(x >tt)], rownames(bio)), ]$group, function(x) unlist(strsplit(x,","))))))
  ccpr <- apply(ccp, 2, function(x) table(unlist(lapply(bio[match(rownames(ccp)[which(x >tt)], rownames(bio)), ]$group, function(x) unlist(strsplit(x,","))))))
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

cg2bed<-function(cg,extend=100){
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

write.bed<-function(bed,file,extend=0){
  bed[,2]<-as.numeric(as.character(bed[,2]))-extend
  bed[,3]<-as.numeric(as.character(bed[,3]))+extend
  if(ncol(bed)==3){
    bed[,4]<-paste(bed[,1],":",bed[,2],"-",bed[,3],sep="")  
  }
  if(ncol(bed)>=4){
    write.table(bed,file=file,sep="\t",col.names=F,row.names=F,quote=F)
  }
}

readmeth450<-function(){
  rlt<-list()
  library("stringr")
  file<-list.files(pattern="jhu*")
  data<-c()
  for(i in file){
    tmp<-read.table(i,head=T,skip=1,row.names=1,sep="\t",check.names = FALSE,as.is=T)
    data<-cbind(data,tmp[,1])
    print(i)
  }
  #load("PancancerMethMatrix_March2016.RData")
  #load("PancancerMethMatrix_March2016.Test.RData")
  colnames(data)<-unlist(lapply(unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[6])),function(x) substr(x,1,15)))
  rownames(data)<-rownames(tmp)
  cancertype<-unique(unlist(lapply(file,function(x) unlist(strsplit(x,"_|.Human"))[2])))
  sampletype<-unlist(lapply(unlist(lapply(file,function(x) unlist(strsplit(x,"[.]"))[6])),function(x) substr(x,14,15)))
  save(data,file=paste(cancertype,"meth.RData",sep="."))
  rlt$data<-data
  rlt$cancertype<-cancertype
  rlt$sampletype<-sampletype
  rlt$cpg<-rownames(data)
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
#gsirlt<-gsirlt[-which(unlist(lapply(gsirlt$group,function(x) length(grep("WBC|Spleen",x))>0))),]
gsirlt1<-subset(gsirlt,refMax>0.6 & GSI>0.65)  # refMax=c(0.3-0.6), GSI=c(0.5-0.7)
gsirlt1$group<-as.character(gsirlt1$group)
bio<-gsi2bio(gsirlt1)

rlt1<-prediction(Data,bio,tt=0.3)              # tissue=c(0.3-0.6)
rlt2<-prediction(data,bio,tt=0.3)              # tissue=c(0.3-0.6)
rlt3<-prediction(rename2(data),bio,tt=0.3)      # reanmed data
sum(rlt1$prediction==names(rlt1$prediction))/length(rlt1$prediction)

# colon, lung and WBC are shared large number common hyermethyalted regions compared with othter tissues.
pdf("heatmap-tovb.pdf")                      # heatmap to all the samples with occurred tissue marker vs background
heatmaprlt<-HeatMap(rlt1$score)
heatmaprlt<-HeatMap(t(rlt1$score))
heatmaprlt<-HeatMap(rlt2$score)
heatmaprlt<-HeatMap(t(r2t1$score))
dev.off()

input<-data.matrix(Data[na.omit(match(rownames(bio),rownames(Data))),])
heatmap(input,cexRow=0.5,cexCol=0.5,output="Performance-of-markers-applied-in-plasma-tissue-of-mapping-in-tissue.pdf")

RRBS<-rrbsdata[,grep("LC-P|CRC-P|NC-P",colnames(rrbsdata))]
input2<-data.matrix(RRBS[na.omit(match(rownames(bio),rownames(RRBS))),])
input2<-RowNARemove(input2,missratio=0.3)
heatmap(input2,cexRow=0.5,cexCol=0.25,output="Performance-of-markers-applied-in-plasma-tissue-of-mapping-in-plasma.pdf")


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
rlt2<-ZscorePredictionPlasma(ccp,ncp,bio,tt=i)    # tt is threshold to binary the MHL matrix, ccp: colon cancer plasma, ncp: normal plasma trainning dataset
prin1<-sum(unlist(apply(rlt1$score,2,function(x) rownames(rlt1$score)[which.max(x)]))=="Colon")
prin2<-sum(unlist(apply(rlt2$score,2,function(x) rownames(rlt1$score)[which.max(x)]))=="Lung")
print(c(i,prin1,prin2))
}




# TCGA dataset, okay, now we can get all the regions overlapped with RRBS, WGBS and 450K




































