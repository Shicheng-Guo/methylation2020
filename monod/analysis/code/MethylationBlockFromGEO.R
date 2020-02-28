#####################################################################################
###   Title : Methylation Block Region analysis based on Methylation 450K beadchip
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Section 1.  function loading
###   Section 2.  Annotation and databased pre-loading 
###   Section 3.  MBR analysis for each GSE dataset
###   Section 4.  Summary 
#####################################################################################
RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
}   

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


table.binomial.test<-function(region,event.time){
  region1<-table(region)
  region2<-unique(region)
  
  p=length(region)/(length(region2)*event.time)
  pvalue<-c()
  for(i in 1:length(table(region))){
    ptmp<-binom.test(table(region)[i],event.time,p,alternative="two.sided")$p.value
    pvalue<-c(pvalue,ptmp)
  }
  
  return(pvalue)
  print(paste("there are",sum(pvalue<0.05/length(pvalue),"regions were over-preferred in the samples"),sep=" "))
}


mbsearch<-function(data,sortedmapfile,window=100,mincpgnumber=4){
  # rowname of data is cpg site and column name is sample id
  # map file is bed file and the fourth column is cpg 
  library("impute")
  output<-list()
  data<-RawNARemove(data)
  data<-impute.knn(data.matrix(data))$data
  #  data<-na.omit(data)
  map<-sortedmapfile
  map<-map[map[,4] %in% rownames(data),]
  newdata<-data[match(map[,4],rownames(data)),]
  if(nrow(map)==nrow(newdata)){
    a<-map[,2]
    i=1
    tmp<-c()
    rlt<-c()
    rowname<-c()
    index<-0
    while(i < length(a)){
      start=a[i]
      end=a[i+1]
      if(end-start<window && end-start>0){
        tmp<-c(tmp,i)
        i=i+1
      }else{
        if(length(tmp)>mincpgnumber){
          index=index+1
          tmp<-c(min(tmp),max(tmp),length(tmp),a[max(tmp)]-a[min(tmp)],round(length(tmp)/(a[max(tmp)]-a[min(tmp)]),4),round((a[max(tmp)]-a[min(tmp)])/(length(tmp))))
          rlt<-rbind(rlt,tmp)
          tmp2<-paste(map[tmp[1],1],":",map[tmp[1],2],"-",map[tmp[2],2],sep="")
          rowname<-c(rowname,tmp2)
        }
        tmp<-c()
        i=i+1
      }
    }
    rownames(rlt)<-rowname
    colnames(rlt)<-c("rowstart","rowstart","CpGNumber","RegionLength","CpGRatio","AverageGap")
    cor<-c()
    for(j in 1:nrow(rlt)){
      cor1<-mean(cor(t(newdata[rlt[j,1]:rlt[j,2],]),use="complete.obs")) # cancer
      cor<-rbind(cor,cor1)
    }
    rownames(cor)<-rowname
    output$hdr<-rlt
    output$hdrc<-cor
    
    return(output)
  }
}


######################################################################################################################
###################################  Annotation and databased pre-loading  ###########################################
######################################################################################################################
library("GEOquery")

setwd("G:\\geo")
map<-read.table("C:/Users/shicheng/Dropbox/Project/methylation/monod/GPL13534.sort.bed",as.is=T,sep="\t") # for windows

setwd("/home/sguo/monod/data/geo")
map<-read.table("/home/sguo/annotation/GPL13534.sort.bed",as.is=T,sep="\t")  # for linux

#######################################################################################################################
##########################  Methylation block region analysis based on GSE36064  ######################################
#####################  78 Healthy Children (peripheral blood leukocytes from healthy children #########################
#######################################################################################################################
library("GEOquery")
load("GSE36064_matrix.Rdata")
data <- as.data.frame(exprs(GSE36064[[1]]))
phen <- pData(phenoData(GSE36064[[1]]))
rlt<-mbsearch(data,sortedmapfile=map)
save(rlt,file="GSE36064.HDR.MB.Rlt.RData")

load("GSE36064.HDR.MB.Rlt.RData")
plot(density(rlt$hdrc),xlim=c(0,1),main="GSE36064")
legend("topright",legend=c("Normal","Schizophrenia"),col=c("red","blue"),lwd=3,bty="n",lty=1,cex=1)

######################################################################################################################
##########################  Methylation block region analysis based on GSE41169  #####################################
######################################################################################################################
library("GEOquery")
load("GSE41169_matrix.Rdata")
data <- as.data.frame(exprs(GSE41169[[1]]))
phen <- pData(phenoData(GSE41169[[1]]))
phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.7),"[:]"),function(x) as.numeric(unlist(x)[2]))  # status 1:control, 2:scz
phen1[phen1==1]<-"Normal"
phen1[phen1==2]<-"Schizophrenia"
data1<-data[,which(phen1=="Normal")]
data2<-data[,which(phen1=="Schizophrenia")]
rlt1<-mbsearch(data1,sortedmapfile=map)
rlt2<-mbsearch(data2,sortedmapfile=map)
rlt<-list()
rlt$normal<-rlt1
rlt$schizophrenia<-rlt2
save(rlt,file="GSE41169.HDR.MB.Rlt.RData")

load("GSE41169.HDR.MB.Rlt.RData")
plot(density(rlt$normal),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$schizophrenia),col="green",lwd=3)
legend("topright",legend=c("Normal","Schizophrenia"),col=c("red","blue"),lwd=3,bty="n",lty=1,cex=1)


######################################################################################################################
##########################  Methylation block region analysis based on GSE42861  #####################################
#############  GSE42861:Differential DNA methylation in 354 Rheumatoid arthritis and 337 normal  #####################
######################################################################################################################
library("GEOquery")
load("GSE42861_matrix.Rdata")
data <- as.data.frame(exprs(GSE42861[[1]]))
phen <- pData(phenoData(GSE42861[[1]]))
phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.1),": "),function(x) unlist(x)[2])  # status
table(phen1)
data1<-data[,which(phen1=="Normal")]
data2<-data[,which(phen1=="rheumatoid arthritis")]
ncol(data1)
ncol(data2)
rlt1<-mbsearch(data1,sortedmapfile=map)
rlt2<-mbsearch(data2,sortedmapfile=map)
rlt<-list()
rlt$normal<-rlt1
rlt$RA<-rlt2
save(rlt,file="GSE42861.HDR.MB.Rlt.RData")

# you can run it if you have saved the result and you can download the data into your PC and plot it in Rstudio
load("GSE42861.HDR.MB.Rlt.RData")
par<-par()
plot(density(rlt$normal,col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$RA),col="green",lwd=3)
legend("topright",legend=c("Normal","Rheumatoid Arthritis"),col=c("red","blue"),lwd=3,bty="n",lty=1,cex=1)

######################################################################################################################
##########################  Methylation block region analysis based on GSE54115  #####################################
######################################################################################################################
library("GEOquery")
load("GSE54115_matrix.Rdata")
data <- as.data.frame(exprs(GSE54115[[1]]))
phen <- pData(phenoData(GSE54115[[1]]))

phen1<-as.character(phen$characteristics_ch1.2)
phen1[1:9]<-"Yamanaka factors"
phen1[10:15]<-"Thomson factors"
phen1[16]<-"ES"
phen1[17:18]<-"Human Fibroblasts"

data1<-data[,which(phen1=="Yamanaka factors")]
data2<-data[,which(phen1=="Thomson factors")]
data3<-data[,which(phen1=="Human Fibroblasts")]

rlt1<-mbsearch(data=data1,sortedmapfile=map)
rlt2<-mbsearch(data=data2,sortedmapfile=map)
rlt3<-mbsearch(data=data3,sortedmapfile=map)

rlt<-list()
rlt$yamanaka<-rlt1
rlt$thomson<-rlt2
rlt$hf<-rlt3
save(rlt,file="GSE54115.HDR.MB.Rlt.RData")

# you can run it if you have saved the result
load("GSE54115.HDR.MB.Rlt.RData")
par<-par()
plot(density(rlt$hf),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$yamanaka),col="green",lwd=3)
lines(density(rlt$thomson),col="blue",lwd=3)
legend("topright",legend=c("Human Fibroblasts","Yamanaka iPS","Thomson iPS"),col=c("red","green","blue"),lwd=3,bty="n",lty=1,cex=1)

######################################################################################################################
##########################  Methylation block region analysis based on GSE54769  #####################################
######################################################################################################################
load("GSE54769_matrix.Rdata")
data <- as.data.frame(exprs(GSE54769[[1]]))
phen <- pData(phenoData(GSE54769[[1]]))

phen1<-as.character(phen$source_name_ch1)
phen1[1:3]<-rep("iPSCs",3)
phen1[4:12]<-rep("iPS-MSCs",9)

rlt1<-mbsearch(data[,1:3],sortedmapfile=map)
rlt2<-mbsearch(data[,4:12],sortedmapfile=map)

rlt<-list()
rlt$ipsc<-rlt1
rlt$ipsMSCs<-rlt2
save(rlt,file="GSE54769.HDR.MB.Rlt.RData")

plot(density(rlt$ipsMSCs),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$ipsc),col="green",lwd=3)
legend("topright",legend=c("iPS-MSCs","iPSCs"),col=c("red","green"),lwd=3,bty="n",lty=1,cex=1)

######################################################################################################################
##########################  Methylation block region analysis based on GSE35069  #####################################
######################################################################################################################
load("GSE35069_matrix.Rdata")
data <- as.data.frame(exprs(GSE35069[[1]]))
phen <- pData(phenoData(GSE35069[[1]]))

phen1<-sapply(strsplit(as.character(phen$title),"_"),function(x) unlist(x)[1])

data1<-data[,which(phen1=="CD14+")]
data2<-data[,which(phen1=="CD19+")]
data3<-data[,which(phen1=="CD4+")]
data4<-data[,which(phen1=="CD56+")]
data5<-data[,which(phen1=="CD8+")]
data6<-data[,which(phen1=="Eos")]
data7<-data[,which(phen1=="Gran")]
data8<-data[,which(phen1=="Neu")]
data9<-data[,which(phen1=="PBMC")]
data10<-data[,which(phen1=="WB")]


rlt1<-mbsearch(data=data1,sortedmapfile=map)
rlt2<-mbsearch(data=data2,sortedmapfile=map)
rlt3<-mbsearch(data=data3,sortedmapfile=map)
rlt4<-mbsearch(data=data4,sortedmapfile=map)
rlt5<-mbsearch(data=data5,sortedmapfile=map)
rlt6<-mbsearch(data=data6,sortedmapfile=map)
rlt7<-mbsearch(data=data7,sortedmapfile=map)
rlt8<-mbsearch(data=data8,sortedmapfile=map)
rlt9<-mbsearch(data=data9,sortedmapfile=map)
rlt10<-mbsearch(data=data10,sortedmapfile=map)

rlt<-list()
rlt$cd14<-rlt1
rlt$cd19<-rlt2
rlt$cd4<-rlt3
rlt$cd56<-rlt4
rlt$cd8<-rlt5
rlt$eos<-rlt6
rlt$gran<-rlt7
rlt$neu<-rlt8
rlt$pbmc<-rlt9
rlt$wb<-rlt10
save(rlt,file="GSE35069.HDR.MB.Rlt.RData")

setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod")
load("GSE35069.HDR.MB.Rlt.RData")
pdf("GSE35069.hdrc.density.pdf")
plot(density(rlt$cd14$hdrc),col=rainbow(10)[3],lty=1,lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$cd19$hdrc),col=rainbow(10)[2],lty=2,lwd=3)
lines(density(rlt$cd4$hdrc),col=rainbow(10)[1],lty=3,lwd=3)
lines(density(rlt$cd56$hdrc),col=rainbow(10)[4],lty=4,lwd=3)
lines(density(rlt$cd8$hdrc),col=rainbow(10)[5],lty=5,lwd=3)
lines(density(rlt$eos$hdrc),col=rainbow(10)[6],lty=6,lwd=3)
lines(density(rlt$gran$hdrc),col=rainbow(10)[7],lty=7,lwd=3)
lines(density(rlt$neu$hdrc),col=rainbow(10)[8],lty=8,lwd=3)
lines(density(rlt$pbmc$hdrc),col=rainbow(10)[9],lty=9,lwd=3)
lines(density(rlt$wb$hdrc),col=rainbow(10)[10],lty=10,lwd=3)
legend("topright",legend=c("CD14+","CD19+","CD4+","CD56+","CD8+","EOS","GRAN","NEU","PBMC","WB"),col=rainbow(10)[c(3,2,1,seq(4,10))],lty=1:10,lwd=3,bty="n",cex=1)
dev.off()


######################################################################################################################
##########################  Summary Analysis to the Methylation Block Regions   #####################################
######################################################################################################################

# Distribution of Average correlation of high CpG density regions in methylation 450K beadchip
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod")
setwd("/home/sguo/monod/data/geo")

jpeg("GSE35069.HDRC.Distribution.Density.jpg")
load("GSE35069.HDR.MB.Rlt.RData")
plot(density(rlt$cd14$hdrc),col=rainbow(10)[3],lty=1,lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$cd19$hdrc),col=rainbow(10)[2],lty=2,lwd=3)
lines(density(rlt$cd4$hdrc),col=rainbow(10)[1],lty=3,lwd=3)
lines(density(rlt$cd56$hdrc),col=rainbow(10)[4],lty=4,lwd=3)
lines(density(rlt$cd8$hdrc),col=rainbow(10)[5],lty=5,lwd=3)
lines(density(rlt$eos$hdrc),col=rainbow(10)[6],lty=6,lwd=3)
lines(density(rlt$gran$hdrc),col=rainbow(10)[7],lty=7,lwd=3)
lines(density(rlt$neu$hdrc),col=rainbow(10)[8],lty=8,lwd=3)
lines(density(rlt$pbmc$hdrc),col=rainbow(10)[9],lty=9,lwd=3)
lines(density(rlt$wb$hdrc),col=rainbow(10)[10],lty=10,lwd=3)
legend("topright",legend=c("CD14+","CD19+","CD4+","CD56+","CD8+","EOS","GRAN","NEU","PBMC","WB"),col=rainbow(10)[c(3,2,1,seq(4,10))],lty=1:10,lwd=3,bty="n",cex=1)
dev.off()

jpeg("GSE42861.HDRC.Distribution.Density.jpg")
load("GSE42861.HDR.MB.Rlt.RData")
plot(density(rlt$normal$hdrc),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$RA$hdrc),col="blue",lwd=3)
legend("topright",legend=c("Normal","RA"),col=c("red","blue"),lwd=3,bty="n",lty=1,cex=1)
dev.off()

jpeg("GSE54769.HDRC.Distribution.Density.jpg")
load("GSE54769.HDR.MB.Rlt.RData")
plot(density(rlt$ipsMSCs$hdrc),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$ipsc$hdrc),col="blue",lwd=3)
legend("topright",legend=c("iPS-MSCs","iPSCs"),col=c("red","blue"),lwd=3,bty="n",lty=1,cex=1)
dev.off()

jpeg("GSE54115.HDRC.Distribution.Density.jpg")
load("GSE54115.HDR.MB.Rlt.RData")
par<-par()
plot(density(rlt$hf),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$yamanaka),col="green",lwd=3)
lines(density(rlt$thomson),col="blue",lwd=3)
legend("topright",legend=c("Human Fibroblasts","Yamanaka iPS","Thomson iPS"),col=c("red","green","blue"),lwd=3,bty="n",lty=1,cex=1)
dev.off()


load("GSE36064.HDR.MB.Rlt.RData")
plot(density(rlt$hdrc),xlim=c(0,1),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
legend("topleft",legend=c("PBMC in Healthy"),col=c("red"),lwd=3,bty="n",lty=1,cex=1)

load("GSE41169.HDR.MB.Rlt.RData")
plot(density(rlt$schizophrenia$hdrc),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$normal$hdrc),col="blue",lwd=3)
legend("topright",legend=c("Schizophrenia","Normal"),col=c("red","blue"),lwd=3,bty="n",lty=1,cex=1)

load("GSE42861.HDR.MB.Rlt.RData")
plot(density(rlt$normal$hdrc),col="red",lwd=3,main="",cex.lab=1.25,cex.axis=1.25)
lines(density(rlt$RA$hdrc),col="blue",lwd=3)
legend("topright",legend=c("Normal","Rheumatoid Arthritis"),col=c("red","blue"),lwd=3,bty="n",lty=1,cex=0.5)

######################################################################################################################
####################################  the Methylation Block Regions   ################################################
######################################################################################################################
load("GSE35069.HDR.MB.Rlt.RData")
r1<-rownames(rlt$cd14$hdrc)[rlt$cd14$hdrc>0.8]
r2<-rownames(rlt$cd19$hdrc)[rlt$cd19$hdrc>0.8]
r3<-rownames(rlt$cd4$hdrc)[rlt$cd4$hdrc>0.8]
r4<-rownames(rlt$cd56$hdrc)[rlt$cd56$hdrc>0.8]
r5<-rownames(rlt$cd8$hdrc)[rlt$cd8$hdrc>0.8]
r6<-rownames(rlt$eos$hdrc)[rlt$eos$hdrc>0.8]
r7<-rownames(rlt$gran$hdrc)[rlt$gran$hdrc>0.8]
r8<-rownames(rlt$neu$hdrc)[rlt$neu$hdrc>0.8]
r9<-rownames(rlt$pbmc$hdrc)[rlt$pbmc$hdrc>0.8]
r10<-rownames(rlt$wb$hdrc)[rlt$wb$hdrc>0.8]
r<-c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
pvalue=table.binomial.test(r,10)
rr1<-cor2bed(r1)
rr2<-cor2bed(r2)
rr3<-cor2bed(r3)
rr4<-cor2bed(r4)
rr5<-cor2bed(r5)
rr6<-cor2bed(r6)
rr7<-cor2bed(r7)
rr8<-cor2bed(r8)
rr9<-cor2bed(r9)
rr10<-cor2bed(r10)
nrow<-c()
rr<-Rbedtools(functionstring="intersectBed",rr1,rr2,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr3,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr4,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr5,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr6,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr7,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr8,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr9,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr10,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
write.table(rr,file="GSE35069.Share.High.HDRC.0.6.bed",sep="\t",col.names=F,row.names=F,quote=F)




load("GSE42861.HDR.MB.Rlt.RData")
r1<-rownames(rlt$normal$hdrc)[rlt$normal$hdrc>0.6]
r2<-rownames(rlt$RA$hdrc)[rlt$RA$hdrc>0.6]
rr1<-cor2bed(r1)
rr2<-cor2bed(r2)
r<-c(r1,r2)
pvalue=table.binomial.test(r,10)
nrow<-c()
rr<-Rbedtools(functionstring="intersectBed",rr1,rr2,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
write.table(rr,file="GSE42861.Share.High.HDRC.0.6.bed",sep="\t",col.names=F,row.names=F,quote=F)




load("GSE54769.HDR.MB.Rlt.RData")
r1<-rownames(rlt$ipsc$hdrc)[rlt$ipsc$hdrc>0.6]
r2<-rownames(rlt$ipsMSCs$hdrc)[rlt$ipsMSCs$hdrc>0.6]
rr1<-cor2bed(r1)
rr2<-cor2bed(r2)
r<-c(r1,r2)
pvalue=table.binomial.test(r,10)
nrow<-c()
rr<-Rbedtools(functionstring="intersectBed",rr1,rr2,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
write.table(rr,file="GSE54769.Share.High.HDRC.0.6.bed",sep="\t",col.names=F,row.names=F,quote=F)




load("GSE54115.HDR.MB.Rlt.RData")
r1<-rownames(rlt$yamanaka$hdrc)[rlt$yamanaka$hdrc>0.6]
r2<-rownames(rlt$thomson$hdrc)[rlt$thomson$hdrc>0.6]
r3<-rownames(rlt$hf$hdrc)[rlt$hf$hdrc>0.6]
rr1<-cor2bed(r1)
rr2<-cor2bed(r2)
rr3<-cor2bed(r3)
r<-c(r1,r2,r3)
pvalue=table.binomial.test(r,3)
min(pvalue)
nrow<-c()
rr<-Rbedtools(functionstring="intersectBed",rr1,rr2,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
rr<-Rbedtools(functionstring="intersectBed",rr,rr3,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
write.table(rr,file="GSE54115.Share.High.HDRC.0.6.bed",sep="\t",col.names=F,row.names=F,quote=F)

dim(rr1)
dim(rr2)
dim(rr3)



load("GSE36064.HDR.MB.Rlt.RData")
r1<-rownames(rlt$hdrc)[rlt$hdrc>0.6]
rr1<-cor2bed(r1)
r<-c(r1)
pvalue=table.binomial.test(r,10)
min(pvalue)
write.table(rr1,file="GSE36064.Share.High.HDRC.0.6.bed",sep="\t",col.names=F,row.names=F,quote=F)




load("GSE41169.HDR.MB.Rlt.RData")
r1<-rownames(rlt$normal$hdrc)[rlt$normal$hdrc>0.6]
r2<-rownames(rlt$schizophrenia$hdrc)[rlt$schizophrenia$hdrc>0.6]
rr1<-cor2bed(r1)
rr2<-cor2bed(r2)
r<-c(r1,r2)
pvalue=table.binomial.test(r,10)
min(pvalue)
nrow<-c()
rr<-Rbedtools(functionstring="intersectBed",rr1,rr2,opt.string="-wa")
nrow<-c(nrow,nrow(rr))
write.table(rr,file="GSE41169.Share.High.HDRC.0.6.bed",sep="\t",col.names=F,row.names=F,quote=F)


######################################################################################################################
####################################  the Methylation Block Regions   ################################################
######################################################################################################################

file<-list.files(pattern="*Share.High.HDRC.0.6.bed.txt")
tmp1<-read.table(file[1],as.is=T)
tmp2<-read.table(file[2],as.is=T)
tmp3<-read.table(file[3],as.is=T)

tmp12<-Rbedtools(functionstring="intersectBed",tmp1,tmp2,opt.string="-wa")
tmp13<-Rbedtools(functionstring="intersectBed",tmp1,tmp3,opt.string="-wa")
tmp23<-Rbedtools(functionstring="intersectBed",tmp2,tmp3,opt.string="-wa")
tmp123<-Rbedtools(functionstring="intersectBed",tmp12,tmp3,opt.string="-wa")

area1=nrow(tmp1)
area2=nrow(tmp2)
area3=nrow(tmp3)

n12=nrow(tmp12)
n13=nrow(tmp13)
n23=nrow(tmp23)
n123=nrow(tmp123)
write.table(tmp123,file="share.hdrc.GSE35069,GSE41169,GSE42861.txt",sep="\t",col.names=F,row.names=F,quote=F)

setwd("/home/sguo/monod/data/geo")
tmp123<-read.table("share.hdrc.GSE35069,GSE41169,GSE42861.txt",sep="\t",as.is=T)
data2<-read.table("conservative.methylation.block.C11.0.6.normal.bed",sep="\t",as.is=T)
data3<-read.table("conservative.methylation.block.C11.0.6.cancer.bed",sep="\t",as.is=T)

rlt1<-Rbedtools(functionstring="intersectBed",tmp123,data2,opt.string="-wa")
rlt2<-Rbedtools(functionstring="intersectBed",tmp123,data3,opt.string="-wa")
rlt3<-Rbedtools(functionstring="intersectBed",data2,data3,opt.string="-wa")
rlt4<-Rbedtools(functionstring="intersectBed",rlt1,data3,opt.string="-wa")

library("VennDiagram")
# Venn for PBMC
area1=138
area2=590
area3=909
n12=123
n13=121
n23=528
n123=119
pdf("a.pdf")
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, cex=3,cat.cex=2,col.lab="white",cat.pos=c(350,0,180),category =c("GSE35069","GSE41169","GSE42861"),col=c(2:4),fill=2:4,lwd=2,ind = T, list.order = 1:3)
dev.off()
? draw.triple.venn
# Venn for PBMC and TCGA
area1=155
area2=119
cross.area=51
draw.pairwise.venn(area1, area2, cross.area, category = c("PBMC","TCGA-Normal"),col=2:3,cat.pos=c(330,0),fill=2:3,cat.cex=2.5,cex=3.5)

area1=155
area2=276
cross.area=55
draw.pairwise.venn(area1, area2, cross.area, category = c("PBMC","TCGA-Cancer"),rotation.degree = 180,col=2:3,cat.pos=c(180,180),fill=2:3,cat.cex=2.5,cex=3.5)

# Venn for PBMC, TCGA-Normal, TCGA-Cancer
area1=155
area2=276
area3=119
n12=113  #nrow(rlt3)
n13=55   #nrow(rlt2)
n23=51   #nrow(rlt1)
n123=43  #nrow(rlt4)
pdf("a.pdf")
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, cex=3,cat.cex=1.5,col.lab="white",cat.pos=c(350,0,180),category =c("TCGA-Normal","TCGA-Cancer","PBMC"),col=c(2:4),fill=2:4,lwd=2,ind = T, list.order = 1:3)
dev.off()





# Venn for iPS 
file<-list.files(pattern="*Share.High.HDRC.0.6.bed.txt")
tmp4<-read.table(file[4],as.is=T)
tmp5<-read.table(file[5],as.is=T)

tmp45<-Rbedtools(functionstring="intersectBed",tmp4,tmp5,opt.string="-wa")


tmp13<-Rbedtools(functionstring="intersectBed",tmp1,tmp3,opt.string="-wa")






