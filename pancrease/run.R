setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/pancrease/")

tmb<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/TMB-1106.csv")
cufflinks<-read.table("cufflinksOut.txt",head=T,sep="\t",check.names = F)
dim(cufflinks)
newcolnames<-unlist(lapply(strsplit(colnames(cufflinks),"FRAL"),function(x) x[1]))
colnames(cufflinks)<-newcolnames
head(tmb)
head(cufflinks)
tail(sort(table(cufflinks$Symbol)))
idv<-unique(unlist(lapply(strsplit(newcolnames,"-"),function(x) x[1])))[2:8]
mean(tmb[,4])
DD<-list()
for(i in 1:length(idv)){
  temp<-cufflinks[,c(grep(idv[i],colnames(cufflinks)))]
  ND<-temp[,grep("-N",colnames(temp))]
  CD<-temp[,-grep("-N",colnames(temp))]
  XX<-CD-ND
  DD[[i]]=XX
}
DD<-data.frame(DD,check.names=F)
head(DD)

tmby<-tmb[match(colnames(DD),tmb$TumoT_Sample_BaTcode),4]

PP<-c()
for(i in 1:nrow(DD)){
  fit=summary(lm(tmby~as.numeric(DD[i,])))
  if(nrow(fit$coefficients)>1){
  P<-fit$coefficients[2,]
  }else{
  P<-c(NA,NA,NA,NA)  
  }
  PP<-rbind(PP,P)
  print(i)
}
head(PP)
out=data.frame(as.character(cufflinks[,1]),PP,FDR=p.adjust(PP[,4],"fdr"))

head(out)
write.table(out,file="TMB-DEG.out.txt",sep="\t",quote = F,col.names = T,row.names = F)


library("qqman")
library("Haplin")
pdf("qqplot.pdf")
pvalues=na.omit(out[,5])
pQQ(pvalues, nlabs =length(pvalues), conf = 0.95)
dev.off()



setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/pancrease/medip")
file=c("2019032901","2019032903","2019040901","2019051703","2019052301","2019053101","2019053102")
FOLD<-list()
for (i in 1:length(file)){
  print(i)
  input<-read.table(paste(file[i],".matrix",sep=""),head=T,row.names=1,check.names=F)
  colnames(input)<-unlist(lapply(strsplit(colnames(input),"_20190"),function(x) x[1]))
  foldc<-data.frame(input[,2:4]/input[,1],check.names = F)
  FOLD[[i]]<-foldc
}
FOLD<-data.frame(FOLD)
foldc<-na.omit(foldc)
foldc<-foldc[-which(!is.finite(foldc))]
xsel<-head(order(foldc,decreasing=T),n=2000)
pdf(paste(i,".hvar.matrix.pdf",sep=""))
temp<-input[match(names(foldc)[xsel],rownames(input)),]
z<-apply(temp,1,function(x) (x-mean(x))/sd(x))
HeatMap(data.matrix(t(z)))
dev.off()


file=list.files(pattern="*.tab")
data<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i])
  data<-cbind(data,temp[,6])
  print(i)
}
colnames(data)<-unlist(lapply(strsplit(file,"_201907"), function(x) x[1]))
rownames(data)<-temp[,1]
file=c("2019032901","2019032903","2019040901","2019051703","2019052301","2019053101","2019053102")
DD<-list()
for (i in 1:length(file)){
  print(i)
  temp<-data[,grep(file[i],colnames(data))]+1
  ND<-temp[,grep("_N",colnames(temp))]
  CD<-temp[,-grep("_N",colnames(temp))]
  head(ND)
  head(CD)
  XX<-CD/ND
  DD[[i]]=XX
}
DD<-data.frame(DD,check.names=F)
tmb2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/TMB_2.txt",head=T)
tmby<-tmb[match(colnames(DD),tmb2$TumoT_Sample_BaTcode),4]
PP<-c()
for(i in 1:nrow(DD)){
  fit=summary(lm(tmby~as.numeric(DD[i,])))
  if(nrow(fit$coefficients)>1){
    P<-fit$coefficients[2,]
  }else{
    P<-c(NA,NA,NA,NA)  
  }
  PP<-rbind(PP,P)
  print(i)
}
head(PP)
out=data.frame(rownames(temp),PP,FDR=p.adjust(PP[,4],"fdr"))
out=out[order(out[,5]),]
head(out)
write.table(out,file="TMB-DMR.out.txt",sep="\t",quote = F,col.names = T,row.names = F)
library("qqman")
library("Haplin")
pdf("qqplot.pdf")
pvalues=na.omit(out[,5])
pQQ(pvalues, nlabs =length(pvalues), conf = 0.95)
dev.off()





HeatMap<-function(data,Rowv=T,Colv=T){
  
  #  note: this function include correlation based heatmap (pearson or spearman)
  #  data: row is gene and column is sample
  #  colname and rowname cannot be NULL  
  #  Usage example:
  #  test<- matrix(runif(100),nrow=20)
  #  colnames(test)=c("A","A","A","B","B")
  #  rownames(test)=paste("Gene",1:20,sep="")
  #  HeatMap(test)
  
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  
  Hclust=function(x){hclust(x,method="complete")}
  Distfun=function(x){as.dist((1 - cor(t(x),method = "spearman")))}
  
  ColSideColors=sidecol(colnames(data))
  
  heatmap.2(data,trace="none", 
            hclust=Hclust,
            distfun=Distfun, 
            cexRow = 1, cexCol = 1,
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=Colv,Rowv = Rowv,
            keysize=0.9, margins = c(5, 10)
  )
}



