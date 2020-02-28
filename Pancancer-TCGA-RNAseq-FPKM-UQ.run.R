##############################################################################################
##############################################################################################
##############################################################################################
manifest="gdc_manifest.2019-03-09.txt"
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste0("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system("curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > file_metadata.txt")

#####################################################################################
############### ENSG, ENST, ENSP, Symbol in hg19 ####################################
#####################################################################################
Symbol1<-read.table("/home/guosa/hpc/temp/Homo_sapiens.GRCh38.95.entrez.tsv",head=T,sep="\t")
Symbol2<-read.table("~/hpc/db/hg19/ENST2Symbol.hg19.txt",head=F,sep="\t")
Symbol2$ENST1<-unlist(lapply(strsplit(as.character(Symbol2$V2),split="[.]"),function(x) x[1]))
Symbol2$ENST2<-unlist(lapply(strsplit(as.character(Symbol2$V2),split="[_]"),function(x) x[1]))
Symbol2<-Symbol2[,c(2:6,13,17:18)]
head(Symbol1)
head(Symbol2)
Symbol<-merge(Symbol1,Symbol2,by.x="transcript_stable_id",by.y="ENST1")
Symbol<-Symbol[,c(11:15,1:3,10,16)]
head(Symbol)
write.table(Symbol,file="ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t",col.names = F,row.names = F,quote=F)
#####################################################################################
id2phen4<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
}

id2phen3<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*"))
}

id2bin<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
  as.numeric(lapply(strsplit(filename,"-"),function(x) x[4]))
}

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}

RawZeroRemove<-function(data,missratio=0.5){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) sum(x==0)>=threshold))
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}
###################################################################################################################
##################### ENSG, ENST, ENSP, Symbol in hg19 ############################################################
###################################################################################################################

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq")
load("rnaseqdata.pancancer.RData")
phen1=read.table("/home/guosa/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("/home/guosa/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
head(phen1)
head(phen2)
head(phen)

colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]

phen1<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen2<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen3<-id2bin(phen$cases.0.samples.0.submitter_id)

phen$bin=phen3

include<-which(c(phen$bin==1 | phen$bin==11))
phen<-phen[include,]
input<-rnaseqdata[,include]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]
colnames(input)<-phen$id
Seq<-paste(phen$project_id,phen$bin,sep="-")

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
library("metafor")
data<-input
data[1:5,1:5]
i=500
Seq<-paste(phen$project_id,phen$bin,sep="-")
mean<-tapply(as.numeric(data[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(data[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(data[i,]),Seq,function(x) length(x))

exclude<-names(which(table(unlist(lapply(strsplit(names(mean),"-"),function(x) x[2])))<2))
exclude <-grep(paste(exclude,collapse="|"),phen$project_id)

phen<-phen[-exclude,]
input<-input[,-exclude]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
colnames(input)<-phen$id
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]

input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
Seq<-paste(phen$project_id,phen$bin,sep="-")

rlt<-c()
coll<-c()
for(i in 1:nrow(input)){
print(i)
mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
m1i=mean[seq(1,length(mean),by=2)]
m2i=mean[seq(2,length(mean),by=2)]
sd1i=sd[seq(1,length(mean),by=2)]
sd2i=sd[seq(2,length(mean),by=2)]
n1i=num[seq(1,length(mean),by=2)]
n2i=num[seq(2,length(mean),by=2)]
Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
output$source=Source
output<-na.omit(output)
es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
coll<-c(coll,i)
}
rownames(rlt)<-rownames(input)[coll]
colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
rlt<-data.frame(rlt)

write.table(rlt,file="TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.txt",sep="\t",quote=F,col.names=NA,row.names=T)
save.image("TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.RData")
################################################################################################################################
################################################################################################################################
################################################################################################################################
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq")
diff<-read.table("TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.txt",head=T,sep="\t")
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
head(diff)
genes <- unlist(lapply(strsplit(as.character(diff[,1]),split="[.]"),function(x) x[1]))
diff$X=genes
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
G_list
output<-merge(diff,G_list,by.x="X",by.y="ensembl_gene_id")
library("Haplin")
pQQ(output$pval, nlabs =nrow(output), conf = 0.95) 
head(output)
write.table(output,file="TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",sep="\t",quote=F,col.names=NA,row.names=T)
################################################################################################################################
##########################  Find Small Proteins  ###############################################################################
################################################################################################################################
library("metafor")
setwd("~/hpc/methylation/Pancancer/RNA-seq")
load("TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.RData")

Symbol1<-read.table("/home/guosa/hpc/temp/Homo_sapiens.GRCh38.95.entrez.tsv",head=T,sep="\t")
Symbol2<-read.table("~/hpc/db/hg19/ENST2Symbol.hg19.txt",head=F,sep="\t")
Symbol2$ENST1<-unlist(lapply(strsplit(as.character(Symbol2$V2),split="[.]"),function(x) x[1]))
Symbol2$ENST2<-unlist(lapply(strsplit(as.character(Symbol2$V2),split="[_]"),function(x) x[1]))
Symbol<-merge(Symbol1,Symbol2,by.x="transcript_stable_id",by.y="ENST1")
NewSymbol<-Symbol[,c(1:4,11:15,22,ncol(Symbol))]
rownames(rlt)<-unlist(lapply(strsplit(as.character(rownames(rlt)),split="[.]"),function(x) x[1]))
Rowname<-NewSymbol[match(rownames(rlt),NewSymbol$gene_stable_id),]$V13


Match<-match(c("ZNF132","TLX2","SOX11"),rownames(input))
CHOL<-c("ADCY2","AHRR","ARRB2","BCAN","DLX5","DMRTA2","FGF8","GRIN1","GSC","HIST1H3G","HOXA1","HOXA9","HOXD12","HOXD4","IFFO1","LHX1","LHX8","LINC00461","MNX1","NR2E1","OTX1","PEG3","PITX2","POU4F3","RBMY1F","SATB2-AS1","SHOX2","SIM1","SLC2A14","SPEG","TFAP2E","TRH","TRIM58","UMOD","ZNF518B")
ESCA2<-c("ACSS3","AKAP12","BCAT1","CD200","CDO1","CYP2A13","DMRTA21","DMRTA22","FEZF2","FNDC1","FOXB1","FOXE1","GHSR","GRIK4","GSC","GSX1","HOXD9","ITGA8","LINC00466","LINE-1","LncRNA","MIR129-2","NKAIN4","NKX1-2","NOTO","NR2E1","NTM","PCDH10","PDE4B","PEX5L","RBMY1F","RUNDC3B","RXFP3","RYR2","SHOX2","SLC6A2","SOX11","SPAG6","SRCIN1","THY1","TLX1","TRH","VSTM2A","WNT7B","ZNF273","ZNF385D","ZNF415","ZNF583","ZNF781")
ESCA1<-c("TRIM71","USP44","SCOC","HMX3","ZNF878","ZNF345","TACC2","ADHFE1","TMEM132C","ZNF385B","GFRA1","DPY19L2P4","cg15830431","KIAA0226L","cg20655070","ZNF71","KCNA6","SLITRK5","KCNA3","LRAT","ELMO1","ZNF793","ZNF542","ZNF844","TLX2","cg05249644","CH25H","ZNF418","ZNF570","TBX4","SALL1","SPATA32","EOMES","ZNF790","ZIK1","LINE-1","cg19396867","APC","chrM","CNR1","ZNF132","EYA4","ZNF461","CHST2","ZNF470","ZNF569","TFPI2","RNLS")
TSG<-c("TP53","DNMT1","DNMT3A","DNMT3B","SUV39H1","FSTL1","HDAC7","STAT3","INHBA","TLR4","CD14")

CHOL[CHOL %in% ESCA1]
CHOL[CHOL %in% ESCA2]

FULL<-c(CHOL,ESCA1,ESCA2,TSG)
Match<-unique(na.omit(match(FULL,rownames(input))))

Match<-na.omit(match(CHOL,rownames(input)))
Match<-na.omit(match(ESCA2,rownames(input)))
Match<-na.omit(na.omit(match(ESCA1,rownames(input))))
Match<-na.omit(na.omit(match(TSG,rownames(input))))

P<-c()
beta<-c()
for(i in Match){
print(rownames(input)[i])
mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
m1i=mean[seq(1,length(mean),by=2)]
m2i=mean[seq(2,length(mean),by=2)]
sd1i=sd[seq(1,length(mean),by=2)]
sd2i=sd[seq(2,length(mean),by=2)]
n1i=num[seq(1,length(mean),by=2)]
n2i=num[seq(2,length(mean),by=2)]
Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
output$source=Source
output<-na.omit(output)
es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
res <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
pdf(paste(rownames(input)[i],".pdf",sep=""))
plot(res)
text(0, -0.1, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (Q = ",
                                            .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                            ", p = ", .(formatC(res$QEp, digits=2, format="E")), "; ", I^2, " = ",
                                            .(formatC(res$I2, digits=1, format="f")), "%)")))

text(0, -0.25, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (beta = ",
                                            .(formatC(res$beta, digits=2, format="f")), ", se = ", .(formatC(res$se, digits=2, format="f")),
                                            ", zval = ", .(formatC(res$zval, digits=2, format="f")), "; ", P, " = ",
                                            .(formatC(res$pval, digits=2, format="E")), ")")))
text(0, -0.4, pos=4, cex=0.75, rownames(input)[i])
dev.off()
P<-c(P,res$pval)
beta<-c(beta,res$beta)
system("mv *.pdf ~/hpc/methylation/PanCancerRNAseq")
}

rlt<-data.frame(beta,P)
rownames(rlt)<-rownames(input)[Match]
newrlt<-rlt[order(rlt$P),]  
write.table(newrlt,file="~/hpc/methylation/PanCancerRNAseq/ESCA1_2_CHOL.DMR.EXP.diff.txt",sep="\t",quote=F,row.names = T,col.names = NA)

GSE23400<-read.table("~/hpc/methylation/esophageal_carcinoma/GEO/GSE23400.txt",head=T,row.names = 1,sep="\t")
GSE23400$B=-GSE23400$B
GSE23400$logFC=-GSE23400$logFC
output<-GSE23400[GSE23400$Gene.symbol %in% FULL,]
write.table(output,file="~/hpc/methylation/PanCancerRNAseq/ESCS.EXP.diff.txt",sep="\t",quote=F,row.names = T,col.names = NA)


####################################################################################################
######################### Volcono Plot ############################################################
####################################################################################################
data1<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lung/TissueLowExpHyperGene.txt",head=T)
rlt<-data.frame(MethQ0=tapply(data1$Q0,data1$GENESYMBOL,mean),MethQ25=tapply(data1$Q25,data1$GENESYMBOL,mean),TPMBlood=tapply(data1$Blood,data1$GENESYMBOL,mean),Blood2TissueRatio=tapply(data1$B2TRatio,data1$GENESYMBOL,mean))
write.table(rlt,file="../../TissueLowExpHyperGene.txt",sep="\t",quote=F,row.names = T,col.names = NA)

data1<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_FPKM_UQ/TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",head=T,sep="\t")
newdata<-data1[order(data1$pval),]
head(newdata$hgnc_symbol,30)

TSG1<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/COSMIC.TSG.hg19.bed",sep="\t")
head(TSG1)
TSG2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/TSGene2.0.txt",sep="\t")
head(TSG2)
rlt1<-data1[data1$hgnc_symbol %in% TSG1$V4,]
rlt2<-data1[data1$hgnc_symbol %in% TSG2$V2,]
head(rlt1)
head(rlt2)

pdf("../../Volcano.TSG.Pancancer.RNAseq.pdf")
par(mfrow=c(2,2))
with(rlt1, plot(beta, -log10(pval), pch=20, main="Volcano plot", xlim=c(-5,5)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(rlt1, pval<.0000005 ), points(beta, -log10(pval), pch=20, col="red"))
with(subset(rlt1, abs(beta)>1), points(beta, -log10(pval), pch=20, col="orange"))
with(subset(rlt1, pval<.0000005 & abs(beta)>1), points(beta, -log10(pval), pch=20, col="green"))
library(calibrate)
with(subset(rlt1, pval<.00005 & abs(beta)>1), textxy(beta, -log10(pval), labs=hgnc_symbol, cex=.5))

with(rlt2, plot(beta, -log10(pval), pch=20, main="Volcano plot", xlim=c(-5,5)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(rlt2, pval<.00000005 ), points(beta, -log10(pval), pch=20, col="red"))
with(subset(rlt2, abs(beta)>2), points(beta, -log10(pval), pch=20, col="orange"))
with(subset(rlt2, pval<.00000005 & abs(beta)>2), points(beta, -log10(pval), pch=20, col="green"))
library(calibrate)
with(subset(rlt2, pval<.00000005 & abs(beta)>2), textxy(beta, -log10(pval), labs=hgnc_symbol, cex=.5))
dev.off()

####################################################################################################
######################### SOX11 ############################################################
####################################################################################################
setwd("~/hpc/methylation/Pancancer/RNA-seq")
load("rnaseqdata.pancancer.RData")

phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
head(phen1)
head(phen2)
head(phen)

colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]

phen1<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen2<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen3<-id2bin(phen$cases.0.samples.0.submitter_id)

phen$bin=phen3

include<-which(c(phen$bin==1 | phen$bin==11))
phen<-phen[include,]
input<-rnaseqdata[,include]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]
colnames(input)<-phen$id
ENST2Symbol<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
ENSG<-unlist(lapply(strsplit(as.character(rownames(input)),split="[.]"),function(x) x[1]))
Symbol<-ENST2Symbol[match(ENSG,ENST2Symbol$V7),5]


SOX11<-input[27814,]
match("SOX11",Symbol)
COR<-c()
for(i in 1:nrow(input)){
  temp<-cor.test(input[27814,],input[i,],use="pairwise.complete.obs")
  COR<-rbind(COR,c(temp$estimate,temp$p.value))
  print(i)
}
rlt<-data.frame(Symbol,cor=COR[,1],pval=COR[,2])
rlt$absCOR<-abs(rlt$cor)
rlt<-na.omit(rlt)
write.table(rlt,file="../../SOX11.correlation.txt",sep="\t",quote=F,row.names = T,col.names = NA)

miR<-rlt[grep("MIR",rlt$Symbol),]
write.table(miR,file="../../SOX11.miRNA.txt",sep="\t",quote=F,row.names = T,col.names = NA)

miRDb<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/esca/SOX11/miRDb.SOX11.txt",head=T,sep="\t")
head(miRDb)

miRDb$miRNA<-lapply(miRDb$miRNA.Name,function(x) gsub("hsa-","",as.character(x)))
miRDb$miRNA<-lapply(miRDb$miRNA,function(x) gsub("-3p|-5p|-","",as.character(x)))
miRDb$miRNA<-lapply(miRDb$miRNA,function(x) gsub("miR","MIR",as.character(x)))
SOX11miRNA<-miR[miR$Symbol %in% miRDb$miRNA,]
write.table(SOX11miRNA,file="../../SOX11.miRNA.miRDB.txt",sep="\t",quote=F,row.names = T,col.names = NA)
MIR221

pdf("../../SOX.miRNA.pdf")
par(mfrow=c(2,2))
plot(log(input[27814,],2),log(input[match("MIR221",Symbol),],2))
plot(log(input[27814,],2),log(input[match("MIR4685",Symbol),],2))
dev.off()
  

tsg<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/TSGene2.0.txt",head=T,sep="\t")
RLT<-subset(rlt,pval<10^-10)
TSG<-subset(RLT[RLT$Symbol %in% tsg$GeneSymbol,],pval<10^-200)
dim(TSG)
sum(TSG$cor<0)
sum(TSG$cor>0)

oncogene<-read.table("ongene_human_803.txt",head=T,sep="\t")
RLT<-subset(rlt,pval<10^-10)
Onco<-subset(RLT[RLT$Symbol %in% oncogene$OncogeneName,],pval<10^-200)
dim(Onco)
sum(Onco$cor<0)
sum(Onco$cor>0)

cor.test(input[27814,],input[match("SALL2",Symbol),],use="pairwise.complete.obs")

###################################################################################################
######################### Survial Analysis to RNA-seq ##########################################
####################################################################################################
library("survival")
library("survminer")

OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
data<-input[,which(id2bin(colnames(input))==1)]
newdata<-data[,na.omit(match(OS$submitter_id,id2phen3(colnames(data))))]
colnames(newdata)<-id2phen3(colnames(newdata))

phen<-OS[match(colnames(newdata),OS$submitter_id),]
head(phen)
phen$censored<-as.numeric(! phen$censored)
phen$week=phen$time/7
head(phen)

ENST2Symbol<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
ENSG<-unlist(lapply(strsplit(as.character(rownames(newdata)),split="[.]"),function(x) x[1]))
Symbol<-ENST2Symbol[match(ENSG,ENST2Symbol$V7),5]

HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
HR1<-HR

HR<-HR1
NewRowName<-paste(Symbol,rownames(newdata),sep="_")
HRV=cbind(HR,Symbol,rownames(newdata))
write.table(HRV,file="../../RNAseq_Suvival.HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)


HRVV<-read.table(file="../../RNAseq_Suvival.HR.txt",sep="\t",head=T)

DATA<-HRVV[which(as.numeric(HRVV[,6])<10^-50),]

for(j in 1:nrow(DATA)){
gene=DATA[j,7]
i=match(gene,Symbol)
dat<-data.frame(Rna=newdata[i,],phen)
dat$Rna[dat$Rna<=median(dat$Rna)]<-0
dat$Rna[dat$Rna>median(dat$Rna)]<-1
coxph(Surv(week,censored)~Rna,dat)
hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
fit <- survfit(Surv(week,censored)~Rna, data = dat)
survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
                  fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                  palette = c("#E7B800","#2E9FDF"),
                  legend = "bottom",legend.title = paste(gene,"HR=",round(hr[2],2),sep=" "),
                  legend.labs = c("Low-Expression","High-Expression"))
ggsave(file = paste("../../Survival_Figure/",gene,".Pancancer.pdf",sep=""), survp$plot)
}

CHOL<-c("ADCY2","AHRR","ARRB2","BCAN","DLX5","DMRTA2","FGF8","GRIN1","GSC","HIST1H3G","HOXA1","HOXA9","HOXD12","HOXD4","IFFO1","LHX1","LHX8","LINC00461","MNX1","NR2E1","OTX1","PEG3","PITX2","POU4F3","RBMY1F","SATB2-AS1","SHOX2","SIM1","SLC2A14","SPEG","TFAP2E","TRH","TRIM58","UMOD","ZNF518B")
ESCA2<-c("ACSS3","AKAP12","BCAT1","CD200","CDO1","CYP2A13","DMRTA21","DMRTA22","FEZF2","FNDC1","FOXB1","FOXE1","GHSR","GRIK4","GSC","GSX1","HOXD9","ITGA8","LINC00466","LINE-1","LncRNA","MIR129-2","NKAIN4","NKX1-2","NOTO","NR2E1","NTM","PCDH10","PDE4B","PEX5L","RBMY1F","RUNDC3B","RXFP3","RYR2","SHOX2","SLC6A2","SOX11","SPAG6","SRCIN1","THY1","TLX1","TRH","VSTM2A","WNT7B","ZNF273","ZNF385D","ZNF415","ZNF583","ZNF781")
ESCA1<-c("TRIM71","USP44","SCOC","HMX3","ZNF878","ZNF345","TACC2","ADHFE1","TMEM132C","ZNF385B","GFRA1","DPY19L2P4","cg15830431","KIAA0226L","cg20655070","ZNF71","KCNA6","SLITRK5","KCNA3","LRAT","ELMO1","ZNF793","ZNF542","ZNF844","TLX2","cg05249644","CH25H","ZNF418","ZNF570","TBX4","SALL1","SPATA32","EOMES","ZNF790","ZIK1","LINE-1","cg19396867","APC","chrM","CNR1","ZNF132","EYA4","ZNF461","CHST2","ZNF470","ZNF569","TFPI2","RNLS")
MFM<-c("ANK1","APC","CDH1","CDH11","CDH13","CLDN11","COL1A2","DAPK","EGFR","ERa","ESR1","FHIT","GATA5","GSTP1","HACE1","hMLH1","HOXA11","HOXA9","KCNH8","L1TD1","LOX","MAPK13","MGMT","CDKN2A ","CDKN2B","PAX6","PTEN","RXRB","RASSF1","RNF","SPAG6","SYK","THBS1","TIMP3","TNFSF10D","ZAR1","ZNF677","RXRA")
GSCG<-c("ERF","DNMT1","DNMT3A","DNMT3B","SUV39H1","FSTL1","INHBA")

FULL<-c(CHOL,ESCA1,ESCA2)
FULL<-c(MFM)
FULL<-c(GSCG)
FULL=FULL[FULL %in% Symbol]
for(j in FULL){
  gene=j
  i=match(gene,Symbol)
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  coxph(Surv(week,censored)~Rna,dat)
  hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
  fit <- survfit(Surv(week,censored)~Rna, data = dat)
  survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
                    fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                    palette = c("#E7B800","#2E9FDF"),
                    legend = "bottom",legend.title = paste(gene,"HR=",round(hr[2],2),sep=" "),
                    legend.labs = c("Low-Expression","High-Expression"))
  ggsave(file = paste("../../Survival_Figure/MFM/",gene,".Pancancer.pdf",sep=""), survp$plot)
}

###############################################################################################################################
################################# Meta-analysis to survival data of RNA-seq ##################################################
##############################################################################################################################
library("survival")
library("survminer")
library("meta")

OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
data<-input[,which(id2bin(colnames(input))==1)]
newdata<-data[,na.omit(match(OS$submitter_id,id2phen3(colnames(data))))]
colnames(newdata)<-id2phen3(colnames(newdata))

phen<-OS[match(colnames(newdata),OS$submitter_id),]
head(phen)
phen$censored<-as.numeric(! phen$censored)
phen$week=phen$time/7
head(phen)

# Meta-analysis of survival data:
logHR <- log(c(0.95, 1.5))
selogHR <- c(0.25, 0.35)
fit<-metagen(logHR, selogHR, sm="HR")
forest(fit)

ENST2Symbol<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
ENSG<-unlist(lapply(strsplit(as.character(rownames(newdata)),split="[.]"),function(x) x[1]))
Symbol<-ENST2Symbol[match(ENSG,ENST2Symbol$V7),5]
###############################################################################################################################
######################### TNM association study to RNA-seq dataset#############################################################
##############################################################################################################################
library("survival")
library("survminer")
library("meta")
setwd("~/hpc/project/LungBrainMetastasis/submission")
gene="ERF"
erf<-read.table("ERF.tsv",head=T,sep="\t")
erf$censored<-abs(as.numeric(erf$censored)-2)
erf$label<-abs(as.numeric(erf$label)-1)
head(erf)
coxph(Surv(time,censored)~label,erf)
hr<-summary(coxph(Surv(time,censored)~label,erf))$coefficients[1,]
fit <- survfit(Surv(time,censored)~label, data = erf)
survp<-ggsurvplot(fit, data = erf,conf.int = F,pval = TRUE,
                  fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                  palette = c("#E7B800","#2E9FDF"),
                  legend = "bottom",legend.title = paste(gene,"HR=",round(hr[2],2),sep=" "),
                  legend.labs = c("Wile","Mutation"))
ggsave(file = paste("./",gene,".Mutation.lungBrainMetastasis.pdf",sep=""), survp$plot)

#######################################################################################################################
############################################# TNM association study to RNA-seq dataset ################################
######################################################################################################################

setwd("/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/GSE112658/RNAseq")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/NatureCommunication/GSE112658/RNAseq")
ENST2SymbolRef<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
head(ENST2SymbolRef)
data<-read.table("GSE112656_Counts_matix.txt",head=T,row.names=1,sep="\t")
head(rownames(data))
ENSG<-unlist(lapply(strsplit(as.character(rownames(data)),split="[.]"),function(x) x[1]))
head(ENSG)
Rownames<-ENST2SymbolRef[match(ENSG,ENST2SymbolRef$V7),5]
newdata=data/colMeans(data)
head(newdata)

P=apply(newdata,1,function(x) t.test(x[1:10],x[11:20],na.rm=T)$p.value)
head(P)
output<-data.frame(Rownames,P)
head(output)
head(rlt,50)

target<-c("FSTL1","FSTL3","IQUB","CDON","FGF10","IL4","TNF","MYD88","IL6")
target<-c("RAB40B","THBS2","LRRC32","GDPD5","XRRA1","SLC20A2","KCNK2","SLC7A11","UVRAG")
target<-c("IL7R","THBS2","LRRC32","GDPD5","XRRA1","SLC20A2","KCNK2","SLC7A11","UVRAG")

FSTL<-data.frame(t(newdata[match(target,output$Rownames),]),type=c(rep("RA",10),rep("OA",10)))
colnames(FSTL)<-c(target,"type")
head(FSTL)
par(mfrow=c(3,3))
for(i in 1:length(target)){
  print(i)
  boxplot(FSTL[,i]~type,FSTL,col=c("red","blue"),ylab=target[i])
}


input<-read.table("candidate.txt",sep="\t",head=T)
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/RA/GSE55457.txt",sep="\t",head=T)
out<-data.frame(input,data[match(input$rightID,data$Gene.symbol),])
x1<-subset(out,regulation.m=="down"  & as.numeric(as.character(logFC))>0)
x2<-subset(out,regulation.m=="up"  & as.numeric(as.character(logFC))<0)
xx<-rbind(x1,x2)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
library("metafor")
library("GEOquery")

GPL96<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/GEO/GPL96.anno.txt",sep="\t",head=T)

GSE55457 <- getGEO("GSE55457")
data1 <- as.data.frame(exprs(GSE55457[[1]]))
phen1 <- pData(phenoData(GSE55457[[1]]))

data <- getGEO("GSE55584")
data2 <- as.data.frame(exprs(data[[1]]))
phen2 <- pData(phenoData(data[[1]]))

data <- getGEO("GSE55235")
data3 <- as.data.frame(exprs(data[[1]]))
phen3 <- pData(phenoData(data[[1]]))

newphen1<-phen1$'clinical status:ch1'
newphen2<-phen2$'clinical status:ch1'
newphen3<-phen3$'disease state:ch1'

newphen1<-gsub("normal control","GSE55457_Normal",newphen1)
newphen1<-gsub("rheumatoid arthritis","GSE55457_Case_RA",newphen1)
newphen1<-gsub("osteoarthritis","GSE55457_Normal",newphen1)

newphen2<-gsub("osteoarthritis","GSE55584_Normal",newphen2)
newphen2<-gsub("rheumatoid arthritis","GSE55584_Case_RA",newphen2)

newphen3<-gsub("healthy control","GSE55235_Normal",newphen3)
newphen3<-gsub("osteoarthritis","GSE55235_Normal",newphen3)
newphen3<-gsub("synovial tissue isolated from osteoarthritic joint","GSE55235_Normal",newphen3)
newphen3<-gsub("rheumatoid arthritis","GSE55235_Case_RA",newphen3)

input<-data.frame(data1,data2,data3)
Seq<-c(newphen1,newphen2,newphen3)

Symbol<-GPL96[match(rownames(input),GPL96$ID),2]

P<-c()
beta<-c()
rlt<-c()
coll<-c()
for(i in 1:nrow(input)){
  print(rownames(input)[i])
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  res <- rma(es,slab=source,method = "REML", measure = "SMD",data=output,verbose=TRUE, digits=5, control=list(maxiter=1000))
  rlt<-rbind(rlt,c(i,-res$beta,res$pval,res$ci.lb,res$ci.ub,res$I2,res$tau2))
  coll<-c(coll,i)
  P<-c(P,res$pval)
  beta<-c(beta,res$beta)
  filename=gsub("/","_",Symbol[i])
  if(res$pval<0.00000000005){
    pdf(paste(filename,".pdf",sep=""))
    plot(res)
    text(0, -0.1, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (Q = ",
                                                .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                                ", p = ", .(formatC(res$QEp, digits=2, format="E")), "; ", I^2, " = ",
                                                .(formatC(res$I2, digits=1, format="f")), "%)")))
    
    text(0, -0.25, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (beta = ",
                                                 .(formatC(res$beta, digits=2, format="f")), ", se = ", .(formatC(res$se, digits=2, format="f")),
                                                 ", zval = ", .(formatC(res$zval, digits=2, format="f")), "; ", P, " = ",
                                                 .(formatC(res$pval, digits=2, format="E")), ")")))
    text(0, -0.4, pos=4, cex=0.75, rownames(input)[i])
    dev.off()
    system("mv *.pdf ./meta")
  }
}

rownames(rlt)<-Symbol
colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
rlt<-data.frame(rlt)
write.table(rlt,file="GSE55457_GSE55584_GSE55235_RA_GPL960_meta.txt",sep="\t",quote=F,col.names=NA,row.names=T)

library("Haplin")
pdf("qqplot.pdf")
pQQ(rlt$pval, nlabs =nrow(output), conf = 0.95) 
dev.off()


################################################################################################################################
############################################  Find Small Proteins  #############################################################
################################################################################################################################
library("metafor")
setwd("~/hpc/methylation/Pancancer/RNA-seq")
load("TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.RData")

Symbol1<-read.table("/home/guosa/hpc/temp/Homo_sapiens.GRCh38.95.entrez.tsv",head=T,sep="\t")
Symbol2<-read.table("~/hpc/db/hg19/ENST2Symbol.hg19.txt",head=F,sep="\t")
proteinLen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/Human_73150_Protein_length.txt",sep="\t")
Symbol2$ENST1<-unlist(lapply(strsplit(as.character(Symbol2$V2),split="[.]"),function(x) x[1]))
Symbol2$ENST2<-unlist(lapply(strsplit(as.character(Symbol2$V2),split="[_]"),function(x) x[1]))
Symbol<-merge(Symbol1,Symbol2,by.x="transcript_stable_id",by.y="ENST1")
NewSymbol<-Symbol[,c(1:4,11:15,22,ncol(Symbol))]
rownames(rlt)<-unlist(lapply(strsplit(as.character(rownames(rlt)),split="[.]"),function(x) x[1]))
newrlt<-subset(rlt,pval<10^-10 & beta<0)
RowNameSymbol<-NewSymbol[match(rownames(newrlt),NewSymbol$gene_stable_id),]$V13
LowExpSmallProtein<-data.frame(newrlt,proteinLen[match(RowNameSymbol,proteinLen[,1]),])
LowExpSmallProtein<-LowExpSmallProtein[order(LowExpSmallProtein$V2),]
write.table(LowExpSmallProtein,file="../../TCGA_Pancancer_LowExpSmallProtein.txt",sep="\t",quote=F,col.names=NA,row.names=T)


FDADrug<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/FDA_approved_drugtarget.txt",head=T,sep="\t")
head(FDADrug)
sort(table(FDADrug$Subcellular.location))









