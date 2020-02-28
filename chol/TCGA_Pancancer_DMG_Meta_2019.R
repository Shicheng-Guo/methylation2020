source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("meta")
library("metafor")
library("survival")
library("survminer")
# i=match("ENSG00000213754.2",rownames(input))
# i=match("ENSG00000231246.1",rownames(input))
# i=grep("ENSG00000280054",rownames(input))
# i=grep("ENSG00000243384",rownames(input))
# i=grep("ENSG00000274370",rownames(input))
# i=grep("ENSG00000266601",rownames(input))
load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")
TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
panc<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/panc.txt",head=T)
phen1=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/File_metadata2.txt",header = T,sep="\t")
head(phen1)
head(phen2)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
# prepare phenotype information
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$phen2<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid<-phen$project_id
head(phen)

idx<-which(phen$phen2==1 | phen$phen2==11)
phen<-phen[idx,]
input<-rnaseqdata[,idx]
idx<-which(phen$pid %in% paste("TCGA-",TCGAProjects,sep=""))
phen<-phen[idx,]
input<-input[,idx]
input[1:5,1:5]
input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
Seq<-paste(phen$project_id,phen$phen2,sep="-")
rlt<-c()
coll<-c()
out1<-c()
preplot<-FALSE
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
  Source<-unlist(lapply(strsplit(names(m1i),"-"),function(x) x[2]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  studlab=unlist(lapply(rownames(output),function(x) unlist(strsplit(x,"-"))[2]))
  coll<-c(coll,i)
  m<-metagen(yi,seTE=vi,data = es,comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="SMD")
  if(preplot && m$pval.random <10^-30){
    print(rownames(input)[i])
    pdf(paste(rownames(input)[i],".SMD.PANC.pdf",sep=""))
    forest(m,leftlabs =studlab,
           lab.e = "Intervention",
           pooled.totals = FALSE,
           smlab = "",studlab=studlab,
           text.random = "Overall effect",
           print.tau2 = FALSE,
           col.diamond = "blue",
           col.diamond.lines = "black",
           col.predict = "red",
           print.I2.ci = TRUE,
           digits.sd = 2,fontsize=8,xlim=c(-6,1))
    dev.off()
  }
  fixedEffect<-c(m$TE.fixed,m$lower.fixed,m$upper.fixed,m$pval.fixed)
  randomEffect<-c(m$TE.random,m$lower.random,m$upper.random,m$pval.random)
  out1<-rbind(out1,c(fixedEffect,randomEffect))
}
colnames(out1)<-c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random","lower.random","upper.random","pval.random")
rownames(out1)<-row.names(input)
out1<-data.frame(out1,ENSG2Symbol(rownames(out1)))
write.table(out1,file="TCGA-SMD-DGE-Meta-Pvalue-2019.txt",sep="\t",col.names = NA,row.names = T,quote=F)
