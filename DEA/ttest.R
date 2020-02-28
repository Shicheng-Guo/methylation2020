setwd("/home/sguo/drug/methyarray");
files1<-list.files(pattern="._analysis.txt")
sample_array<-array(substr(files1,34,45))
nrows=27578
ncols=233

cpgname<-(read.table(files1[1],head=T,sep="\t")[,2])
genename<-read.table(files1[1],head=T,sep="\t")[,4]
gene<-matrix(rep(NA,nrows*ncols),nrows,ncols)
for (j in 1:233){
gene[,j]<-unlist(read.table(files1[j],head=T,sep="\t")[,3])
}

load("sample233.RData")
type<-array(data.matrix(data.frame(sample1[match(sample_array,sample1[,1]),13]))-1)

res<-which(type==1)
sen<-which(type==2)

narow<-apply(gene,1,function(x) any(all(is.na(x[res]))), all(is.na(x[sen])))
narows<-which(narow==T)
gene<-gene[-narows,]
cpgname<-cpgname[-narows]
genename<-genename[-narows]
result<-list()
result<-apply(gene,1,function(x) t.test(x[res],x[sen]))

pvalue=array()
effect=array()
for (i in 1:dim(gene)[1]){
pvalue[i]<-result[[i]]$p.value
effect[i]<-result[[i]]$estimate[1]/result[[i]]$estimate[2]  #res increase if great than zero
}

pvalue.fdr<-p.adjust(pvalue, method = "fdr", n = length(pvalue))
pvalue.bonferroni<-p.adjust(pvalue, method = "bonferroni", n = length(pvalue))

output<-data.frame(cpgname,genename,pvalue,pvalue.fdr,pvalue.bonferroni,effect)

output2<-output[order(output[,3]),]
write.table(output2,file="output.txt",sep="\t",col.names=NA,row.names=T,quote=F)

plosgene<-read.table("drug_mnotu_feature.bed",sep="\t",as.is=F)[,4]
drugexp<-output2[match(plosgene,output2[,1]),]
write.table(drugexp,file="drug_expression.txt",col.names=NA,row.names=T,sep="\t",quote=F)

heatmapdata2<-t(gene[match(array(output2[1:15,1]),genename),])
group<-type
gheatmap(heatmapdata2,group,17)






