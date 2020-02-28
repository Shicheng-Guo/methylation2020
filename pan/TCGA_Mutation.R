BiocManager::install("TCGAbiolinks")
BiocManager::install("DT")

library("TCGAbiolinks")
library("DT")
library("maftools")
library("dplyr")

pid<-TCGAbiolinks:::getGDCprojects()$project_id
pid<-pid[grep("TCGA",pid)]
pid<-unlist(lapply(strsplit(pid,"-"),function(x) x[2]))

for(i in pid){
query <- GDCquery_Maf(i, pipelines = "muse")
newdata<-data.frame(query$Tumor_Sample_Barcode,query$Hugo_Symbol,query$Variant_Classification,query$Variant_Type)
head(newdata)
write.table(newdata,file=paste("TCGA-mutation.",i,".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
}

file<-list.files(pattern="TCGA*")
ID<-c()
for(i in file){
  id<-read.table(i,head=T,sep="\t")
  ID<-c(ID,as.character(id[,1]))
}
fid<-read.table("id.txt",sep="\t",head=T)
fid<-fid[fid[,1] %in% id2phen3(ID),]

write.table(fid,file="id.1187.txt",sep="\t",col.names=T,row.names = F,quote=F)

maf <- GDCquery_Maf("CHOL", pipelines = "muse")
maf[1:20,]
datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           legacy = TRUE)


maf <- GDCquery_Maf("CHOL", pipelines = "muse") %>% read.maf
datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = titv)



source("GscTools.R")
id1187<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AIDrugResponse/master/extdata/overlapid/id.1187.txt",head=T,sep="\t")
image<-read.table("gdc_manifest.2019-09-12_image.txt",head=T)
id1109<-id1187[id1187[,1] %in% id2phen3(image$filename),]
write.table(id1109,file="id1109.txt",sep="\t",quote=F,col.names = T,row.names = F)
