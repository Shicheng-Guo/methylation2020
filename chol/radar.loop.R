library("fmsb")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-1")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
data[1:12,1:12]
methdata<-data.matrix(data[,c(7:180)])
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rownames(methdata))
head(colnames(methdata))
methdata[1:5,1:5]

for( Symbol in names(table(unlist(lapply(strsplit(rownames(methdata),"_"),function(x) x[1]))))){
print(Symbol)
input<-methdata[grep(Symbol,rownames(methdata)),]
meth<-apply(input,1,function(x) tapply(x,as.numeric(colnames(input))%%2,function(y) mean(y,na.rm=T)))
meth.sd<-apply(input,1,function(x) tapply(x,as.numeric(colnames(input))%%2,function(y) sd(y,na.rm=T)))
colnames(meth)<-unlist(lapply(strsplit(colnames(meth),"_"),function(x) x[2]))
rd <-rbind(rep(max(meth)+0.1,ncol(meth)) , rep(0,ncol(meth)),meth)
rownames(rd)<-1:4
rd<-data.frame(rd,check.names = F)
head(rd)
colors_border=c(rgb(0.2,0.5,0.5,0.9),rgb(0.8,0.2,0.5,0.9))
colors_in=c(rgb(0.2,0.5,0.5,0.4),rgb(0.8,0.2,0.5,0.4))
pdf(paste("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-1/radar/",Symbol,".radar.pdf",sep=""))
radarchart(rd, axistype=1, 
            pcol=colors_border,pfcol=colors_in,plwd=4,plty=1 , 
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,max(meth)+0.1,0.1), cglwd=1.1,
            vlcex=0.8 
)
legend(x=1, y=1, legend = c("Normal", "Cancer"), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1, pt.cex=1.6)
dev.off()
}



setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-2")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))






