## loading packages
# BiocManager::install("ChIPseeker")
# BiocManager::install("clusterProfiler")
# BiocManager::install("DO.db")
# BiocManager::install("ReactomePA")
BiocManager::install("org.Hs.eg.db")

library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(org.Hs.eg.db)



library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/pancrease/medip")

files <- list.files(pattern="2019032901_*")
peak <- readPeakFile(files[[1]],head=F)
covplot(peak, weightCol="V5")
head(peak)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagMatrix <- tagMatrixList[[4]]
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
peakHeatmap(files[[4]], TxDb=txdb, upstream=3000, downstream=3000, color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf2(files[[4]], TxDb=txdb, upstream=3000, downstream=3000,xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
peakAnno <- annotatePeak(files[[4]], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding loci\nrelative to TSS")
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,fun= "enrichKEGG",pvalueCutoff = 0.05,pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)
p <- GRanges(seqnames=c("chr1", "chr3"),ranges=IRanges(start=c(1, 100), end=c(50, 130)))
shuffle(p, TxDb=txdb)
enrichPeakOverlap(queryPeak= files[[5]],targetPeak=unlist(files[1:4]),TxDb=txdb,pAdjustMethod ="BH",nShuffle= 50,chainFile= NULL,verbose= FALSE)


files <- list.files(pattern="*narrowPeak")
promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)
id<-unique(unlist(lapply(strsplit(basename(files),"_"),function(x) x[1])))
for(i in 1:length(id)){
  file=list.files(pattern=paste(id[i],"_",sep=""))
  tagMatrixList <- lapply(file, getTagMatrix, windows=promoter)
  names(tagMatrixList) <- unlist(lapply(strsplit(basename(file),"_2019"),function(x) x[1]))
  pp<-plotAvgProf(tagMatrixList, xlim=c(-5000, 5000))
  ggsave(paste(id[i],".tss.pdf",sep=""))

  pdf(paste(id[i],".feature.pdf",sep=""))
  peakAnnoList <- lapply(file, annotatePeak,TxDb=txdb,tssRegion=c(-5000, 5000), verbose=FALSE)
  names(peakAnnoList) <- unlist(lapply(strsplit(basename(file),"_2019"),function(x) x[1]))
  plotAnnoBar(peakAnnoList)
  dev.off()
  print(i)
}

for(i in 1:length(id)){
  file=list.files(pattern=paste(id[i],"_",sep=""))
  pdf(paste(id[i],".feature.pdf",sep=""))
  peakAnnoList <- lapply(file, annotatePeak,TxDb=txdb,tssRegion=c(-5000, 5000), verbose=FALSE)
  names(peakAnnoList) <- unlist(lapply(strsplit(basename(file),"_2019"),function(x) x[1]))
  plotAnnoBar(peakAnnoList)
  ggsave(paste(id[i],".feature.pdf",sep=""))
  print(i)
}


#### summary peaks location
files <- list.files(pattern="*narrowPeak")
peak<-matrix(nrow=13,ncol=length(files))
for(i in 1:length(files)){
peakAnno <- annotatePeak(files[i], tssRegion=c(-5000, 5000),TxDb=txdb, annoDb="org.Hs.eg.db")
peak[,i]<-peakAnno@annoStat[,2]
}
colnames(peak)<-unlist(lapply(strsplit(files,"_2019"),function(x) x[1]))
rownames(peak)<-peakAnno@annoStat[,1]
write.table(peak,file="peakDistribution.5k.txt",sep="\t",col.names=NA,row.names=T,quote = F)






#### summary peaks location
pdf(file=paste(files[i],".sinpeak.pdf",sep=""))
vennpie(peakAnno)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
upsetplot(peakAnno, vennpie=F)
dev.off()


