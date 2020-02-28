

library(gplots)
pdf("hap_heatmap.pdf", height=6, width=4)
x=read.table("tmp_hapMatrix.txt",header=TRUE,row.names=1);
heatmap.2(as.matrix(x), col=colorRampPalette(c("steelblue1","purple"))(16), scale="none", cexCol=1.0, cexRow=0.2, trace="row", tracecol="grey", dendrogram="none", Rowv=FALSE,Colv=FALSE,key=F) 
dev.off()

