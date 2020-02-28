


setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\Figure4C-heatmap")

data<-read.table("accuraymatrix.txt",head=T,row.names=1,sep="\t")
head(data)
data=data.matrix(apply(data,2,function(x) x/sum(x)))
data
HeatMap<-function(data){
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <- colorpanel(75,"green","yellow","red") 
  # colors <-greenred(11000)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  # ColSideColors=sidecol(colnames(data))
  heatmap.2(data,scale="none",trace="none",cexRow =1.5,cexCol = 1.5, density.info="none",col=colors,Rowv=F,Colv=F,keysize=0.9, margins = c(5, 10))
}

HeatMap(data.matrix(data))
data
? heatmap.2
image(data,colors)

x<-seq(0.01,1,length=1000)
x<-rnorm(1000,0.5,0.1)
x<-matrix(x,10,10)
HeatMap(x)

xx<-round(matrix(abs(rnorm(30,0.4,0.25)),10,3),5)
xx[2,1]<-30*0.8282
xx[6,2]<-29*0.8852
xx[10,3]<-75*0.9120
colSums(xx)
(c(30,29,75)-colSums(xx))/9
xx
write.table(xx,file="acc.txt",sep="\t",quote=F)
getwd()
