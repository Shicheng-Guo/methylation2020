# GPL13534	Illumina HumanMethylation450 BeadChip (HumanMethylation450_15017482)

library("GEOquery")
destdir="/gpfs/home/guosa/hpc/methylation/GEO"
GSE41169 <- getGEO("GSE41169")
GSE53045 <- getGEO("GSE53045")
GSE35069 <- getGEO("GSE35069")
GSE32148 <- getGEO("GSE32148")
GSE36054 <- getGEO("GSE36054")
GSE36064 <- getGEO("GSE36064")
GSE42861 <- getGEO("GSE42861")
GSE71841 <- getGEO("GSE71841")
GSE35069 <- getGEO("GSE35069")
GSE87571 <- getGEO("GSE87571")

system("wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76269/suppl/GSE76269_betas.extra.epigenomics.tsv.gz")
485513
485577

data1 <- as.data.frame(exprs(GSE41169[[1]]))
data2 <- as.data.frame(exprs(GSE53045[[1]]))
data3 <- as.data.frame(exprs(GSE35069[[1]]))
data4 <- as.data.frame(exprs(GSE32148[[1]]))
data5 <- as.data.frame(exprs(GSE36054[[1]]))
data6 <- as.data.frame(exprs(GSE36064[[1]]))
data7 <- as.data.frame(exprs(GSE42861[[1]]))
data9 <- as.data.frame(exprs(GSE71841[[1]]))
data10 <- as.data.frame(exprs(GSE35069[[1]]))

head(data1[1:5,1:5])
head(data2[1:5,1:5])
head(data3[1:5,1:5])
head(data4[1:5,1:5])
head(data5[1:5,1:5])
head(data6[1:5,1:5])
head(data7[1:5,1:5])
head(data8[1:5,1:5])
head(data9[1:5,1:5])
head(data10[1:5,1:5])

dim(data1)
dim(data2)
dim(data3)
dim(data4)
dim(data5)
dim(data6)
dim(data7)
dim(data8)
dim(data9)
dim(data10)

data<-cbind(data1,data2,data3,data4,data5,data6,data7,data9,data10)
tmp<-t(apply(data,1,function(x) c(cpg=rownames(x),mean=mean(x,na.rm=T),
                                  median=median(x,na.rm=T),
                                  SD=sd(x,na.rm=T),
                                  Quantile=quantile(x,na.rm=T),
                                  SampleSize=length(na.omit(x)))))
write.table(tmp,file="Normal.PBMC.GEO.HM450K.Beta.txt",col.names=NA,row.names=T,sep="\t",quote=F)

normalpbmc450beta<-tmp
save(normalpbmc450beta,file="Normal.PBMC.GEO.HM450K.Beta.RData")

load("/mnt/bigdata/Genetic/Projects/shg047/methylation/GEO/Normal.PBMC.GEO.HM450K.Beta.RData")
normalpbmc450beta

tmp<-data.frame(tmp)
hypo<-subset(tmp,mean<0.3 & median<0.3 & Quantile.75.<0.3)
hype<-subset(tmp,mean>0.6 & median>0.6 & Quantile.75.>0.6)

485513
485577




system("wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76269/suppl/GSE76269_betas.extra.epigenomics.tsv.gz")
system("gunzip GSE76269_betas.extra.epigenomics.tsv.gz")
data<-read.table("~/hpc/methylation/GEO/GSE76269_betas.extra.epigenomics.tsv",head=T,row.names=1,sep="\t")

library("GEOquery")
GSE76269 <- getGEO("GSE76269",getGPL = F)



