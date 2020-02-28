setwd("/home/gsc/tcga/ExpressionGenes/AgilentG4502A_07_3");
#load total array data and then extract the data of patients with explicit drug information
# core 233 samples
#select sample with both bam file and drug response information
#Level: drug info-> bam file
load("fpca.result.RData")
bam_sample<-unique(array(substr(rownames(fpca.rnaseq.rlt),1,12)))   #get the sample order of the cluster analysis
which(bam_sample=="TCGA-23-1023")
nature_sample<-read.table("clinical.txt",head=T,as.is=F,sep="\t")
select<-nature_sample[match(bam_sample,nature_sample$BCRPATIENTBARCODE),]
selected.res<-select[which(select[,13]=="Resistant"),]
selected.sen<-select[which(select[,13]=="Sensitive"),]
dim(selected.res)
dim(selected.sen)
sample1<-rbind(selected.res,selected.sen)
save(sample1,file="sample233.RData")
#Level: drug info-> bam file -> array_expression
files<-list.files(pattern=".txt")
array_samples<-array(substr(files,30,41))
length(array(na.omit(match(unique(array_samples),sample1[,1]))))
#level methylation
setwd("/home/gsc/tcga/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3")
methy_files<-list.files(pattern=".txt")
methy_samples<-array(substr(methy_files,34,45))
methy_samples
length(array(na.omit(match(unique(methy_samples),sample1[,1]))))
#level cnv

##
write.table(sample1,file="core_233_sample.txt",sep="\t",quote=F,)

setwd("/home/gsc/tcga/f0151d84-2235-4e4b-8a78-ef7ed873b21b/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3")
files1<-list.files(pattern=".txt")
array_samples1<-array(substr(files1,30,41))
length(array(na.omit(match(unique(array_samples1),sample1[,1]))))

setwd("/home/gsc/tcga/f0151d84-2235-4e4b-8a78-ef7ed873b21b/Expression-Genes/UNC__AgilentG4502A_07_2/Level_3")
files2<-list.files(pattern=".txt")
array_samples2<-array(substr(files2,30,41))
length(array(na.omit(match(unique(array_samples2),sample1[,1]))))

match(array_samples2,array_samples1)



bam_sample<-(array(substr(rownames(fpca.rnaseq.rlt),1,12)))   #get the sample order of the cluster analysis
match(nature_sample$BCRPATIENTBARCODE,bam_sample)
files<-list.files(pattern=".txt")

load("fpca.result.RData")
x1<-substr(rownames(fpca.rnaseq.rlt),1,16)
x2<-substr(files,30,45)

x<-match(unlist(x1),unlist(x2))
sum(is.na(x))
files2<-files[x]
x3<-substr(files2,30,45)     # all the samples

data<-read.table("unc.edu__AgilentG4502A_07_3__TCGA-61-2096-01A-01R-0670-07__gene_expression_analysis.txt",head=T,sep="\t")
data[,2]
genename<-read.table("drug_asso_gene.txt",head=F,as.is=F,sep="\t")
head(genename)
row<-match(unlist(genename),unlist(data[,2]))[which(is.na(match(unlist(genename),unlist(data[,2])))==F)]
length(row)
length(files)
gene<-matrix(NA,length(row),length(x))
for (j in 1:length(x)){
  gene[,j]<-read.table(files2[j],head=T,sep="\t",as.is=F)[row,3]
  print (j);
}
dim(gene)
sample<-substr(x3,1,12)
colnames(gene)<-sample
rownames(gene)<-unlist(genename)
write.table(gene,file="arraydata.txt",col.names=NA,row.names=T,sep="\t",quote=F)

#check the sample of bam-analysis
#use 332 sample
source("fclust.R")
source("heatmap2.R")
library("fpc")
sampleinfo<-read.table("bam_follow_patient_sample.txt",sep="\t",head=T)
x1<-sampleinfo$Patients
load("fpca.result.RData")
x2<-array(substr(rownames(fpca.rnaseq.rlt),1,12))    #get the sample order of the cluster analysis
group0<-array(data.matrix(data.frame(sampleinfo[match(x2,x1),6])))  # get all the patients info according to the order of cluster analysis
group1<-group0[-which(is.na(group0))]   # remove the patients who info is "NA"
xx332<-x1[match(x2,x1)[-which(is.na(sampleinfo[match(x2,x1),6])==T)]]

##sample
xx332
sum(is.na(match(xx332,sample)))  # sample of bam fpca are all locate in array data
x4<-sample[match(xx332,sample)]
dim(gene)
dim(ll)
group0<-array(data.matrix(data.frame(sampleinfo[match(x4,x1),6])))  # get all the patients info according to the order of cluster analysis
gene2<-gene[,match(xx332,sample)]
gene3<-t(gene2)
dim(gene3)
gene3[which(gene3>13)]<-NA

gheatmap(gene3,rowname=group0,17)

fhclust(gene3,method="euclidean",cut=4,bootstrap=FALSE)
d <- dist(gene3, method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
# plot(fit) # display dendogram
clusters <- cutree(fit, k=4) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 cluster
# rect.hclust(fit, k=cut, border="red") 
list<-names(table(clusters))
rlt<-list()
rlt<-lapply(list,function(x) names(which(clusters==x)))
rlt


d<-data.frame(sample=names(clusters),clusters)
write.table(d,file="result.txt",col.names=T,row.names=F,sep="\t")
d

#succus

#mappedreads<-data2[match(substr(rownames(fpca.rnaseq.rlt),1,51),substr(data2[,1],1,51)),2]
#max<-max(data2[match(substr(rownames(fpca.rnaseq.rlt),1,51),substr(data2[,1],1,51)),2])
#fator<-mappedreads/max
#adjust.factor<-data.frame(substr(rownames(fpca.rnaseq.rlt),1,28),fator)




