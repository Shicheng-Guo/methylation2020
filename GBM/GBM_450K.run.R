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


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/GBM/GEO")
library("GEOquery")

GSE41826 <- getGEO("GSE41826")
data <- as.data.frame(exprs(GSE41826[[1]]))
phen <- pData(phenoData(GSE41826[[1]]))
GSE41826<-list()
GSE41826$beta=data
GSE41826$phen=phen
save(GSE41826,file="GSE41826.RData")

GSE89707 <- getGEO("GSE89707")
data <- as.data.frame(exprs(GSE89707[[1]]))
phen <- pData(phenoData(GSE89707[[1]]))
GSE89707<-list()
GSE89707$beta=data
GSE89707$phen=phen
save(GSE89707,file="GSE89707.RData")

GSE74486 <- getGEO("GSE74486")
data <- as.data.frame(exprs(GSE74486[[1]]))
phen <- pData(phenoData(GSE74486[[1]]))
GSE74486<-list()
GSE74486$beta=data
GSE74486$phen=phen
save(GSE74486,file="GSE74486.RData")

GSE103659 <- getGEO("GSE103659")
data <- as.data.frame(exprs(GSE103659[[1]]))
phen <- pData(phenoData(GSE103659[[1]]))
GSE103659<-list()
GSE103659$beta=data
GSE103659$phen=phen
save(GSE103659,file="GSE103659.RData")

GSE114534 <- getGEO("GSE114534")
data <- as.data.frame(exprs(GSE114534[[1]]))
phen <- pData(phenoData(GSE114534[[1]]))
GSE114534<-list()
GSE114534$beta=data
GSE114534$phen=phen
save(GSE114534,file="GSE114534.RData")

GSE66351 <- getGEO("GSE66351")
data <- as.data.frame(exprs(GSE66351[[1]]))
phen <- pData(phenoData(GSE66351[[1]]))
GSE66351<-list()
GSE66351$beta=data
GSE66351$phen=phen
save(GSE66351,file="GSE66351.RData")


dim(GSE41826$beta)
dim(GSE89707$beta)
dim(GSE74486$beta)
dim(GSE103659$beta)
dim(GSE114534$beta)
dim(GSE66351$beta)

dim(GSE41826$phen)
dim(GSE89707$phen)
dim(GSE74486$phen)
dim(GSE103659$phen)
dim(GSE114534$phen)
dim(GSE66351$phen)

load("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/methdata.pancancer.RData")
methdata[1:5,1:5]
phen4<-id2phen4(colnames(methdata))
phen3<-id2phen3(colnames(methdata))
bin<-id2bin(colnames(methdata))
pid<-id2pid(colnames(methdata))
phen<-data.frame(phen4=phen4,phen3=phen3,pid=pid,bin=bin)
exclude<-which(c(phen$bin !=1 & phen$bin !=11))
phen<-phen[-exclude,]
input<-methdata[,-exclude]
Seq<-paste(phen$pid,phen$bin,sep="-")
head(phen)
input[1:5,1:5]
Seq<-paste(phen$pid,phen$bin,sep="-")
head(phen)
input[1:5,1:5]
BC<-grep("LGG|GBM",colnames(input))
newinput<-input[,BC]
newphen<-phen[BC,]
newinput[1:5,1:5]
head(newphen)

clinical<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/clinical.tsv",head=T,sep="\t")
clinical$age=2019-as.numeric(as.character(clinical$year_of_birth))
clinical$gender=clinical$gender
clinical$tissue=clinical$tissue_or_organ_of_origin
phen4<-id2phen4(colnames(newinput))
phen3<-id2phen3(colnames(newinput))
bid<-id2bin(colnames(newinput))
pid<-id2pid(colnames(newinput))
age<-clinical[match(phen3,clinical$submitter_id),]$age
gender<-clinical[match(phen3,clinical$submitter_id),]$gender
disease=paste(pid,bid,sep="_")
phen<-data.frame(age,disease,gender,tissue="Brain")
rownames(phen)<-phen4
write.table(phen,file="TCGA.clinical.phenotype.txt",sep="\t",col.names=NA,row.names=T,quote=F)

load("GSE41826.RData")
load("GSE89707.RData")
load("GSE74486.RData")
load("GSE103659.RData")
load("GSE114534.RData")
load("GSE66351.RData")
GSE41826$beta<-GSE41826$beta[match(rownames(newinput),rownames(GSE41826$beta)),]
GSE89707$beta<-GSE89707$beta[match(rownames(newinput),rownames(GSE89707$beta)),]
GSE74486$beta<-GSE74486$beta[match(rownames(newinput),rownames(GSE74486$beta)),]
GSE103659$beta<-GSE103659$beta[match(rownames(newinput),rownames(GSE103659$beta)),]
GSE114534$beta<-GSE114534$beta[match(rownames(newinput),rownames(GSE114534$beta)),]
GSE66351$beta<-GSE66351$beta[match(rownames(newinput),rownames(GSE66351$beta)),]
dim(GSE41826$beta)
dim(GSE89707$beta)
dim(GSE74486$beta)
dim(GSE103659$beta)
dim(GSE114534$beta)
dim(GSE66351$beta)
GSE41826$phen<-data.frame(GSE41826$phen[,c(36,38,41)],tissue="brain",dataset="GSE41826")
GSE89707$phen<-data.frame(GSE89707$phen[,c(35,36,37)],tissue="brain",dataset="GSE89707")
GSE74486$phen<-data.frame(age="NA",GSE74486$phen[,c(45,46,47)],dataset="GSE74486")
GSE103659$phen<-data.frame(GSE103659$phen[,c(1,40)],disease="GBM",gender="NA",tissue="brain",dataset="GSE103659")
GSE103659$phen<-GSE103659$phen[,2:6]
GSE114534$phen<-data.frame(age="NA",GSE114534$phen[,c(8,34)],tissue="brain",dataset="GSE114534")
GSE66351$phen<-data.frame(GSE66351$phen[,c(41,45,49)],tissue="brain",dataset="GSE66351")
colnames(GSE41826$phen)=c("age","disease","gender","tissue","dataset")
colnames(GSE89707$phen)=c("age","disease","gender","tissue","dataset")
colnames(GSE74486$phen)=c("age","disease","gender","tissue","dataset")
colnames(GSE103659$phen)=c("age","disease","gender","tissue","dataset")
colnames(GSE114534$phen)=c("age","disease","gender","tissue","dataset")
colnames(GSE66351$phen)=c("age","disease","gender","tissue","dataset")

rownames(GSE41826$phen)
rownames(GSE89707$phen)
rownames(GSE74486$phen)
rownames(GSE103659$phen)
rownames(GSE114534$phen)
rownames(GSE66351$phen)


phen1<-rbind(GSE41826$phen,GSE89707$phen,GSE74486$phen,GSE103659$phen,GSE114534$phen,GSE66351$phen)
phen2<-read.table(file="TCGA.clinical.phenotype.txt",sep="\t",head=T,row.names = 1)
phen2$dataset="TCGA"
phen<-rbind(phen1,phen2)

phen$gender[phen$gender=="female"]<-1
phen$gender[phen$gender=="F"]<-1
phen$gender[phen$gender=="f"]<-1
phen$gender[phen$gender=="Female"]<-1
phen$gender[phen$gender=="Male"]<-0
phen$gender[phen$gender=="M"]<-0
phen$gender[phen$gender=="m"]<-0
phen$gender[phen$gender=="male"]<-0
phen$gender[phen$gender==""]<-NA
phen$gender[phen$gender=="--"]<-NA

phen$tissue[phen$tissue=="Brain"]<-"brain"
phen$tissue[phen$tissue=="Cerebellar cortex"]<-"brain"
phen$tissue[phen$tissue=="Cerebrum"]<-"brain"
phen$tissue[phen$tissue=="Frontal Cortex"]<-"brain"

phen$disease[grep("Astrocytoma|astrocytoma|LGG._1|Ganglioneuroblastom|Hirnstammgliom|glialer",phen$disease)]<-"LGG"
phen$disease[grep("schizophrenia|control|Control|GBM._11|CTRL|Depression|AD|Alzheimer|syndrome",phen$disease)]<-"Normal"
phen$disease[grep("Glioblastoma|GBM._1|Gliosarcoma",phen$disease)]<-"GBM"

table(phen$age)
table(phen$disease)
table(phen$gender)
table(phen$tissue)

write.table(phen,file="MCRI.Brain.Project.phenotype.txt",sep="\t",col.names=NA,row.names=T,quote=F)
phen<-read.table(file="MCRI.Brain.Project.phenotype.txt",head=T,row.names = 1)

newphen<-phen
colnames(newinput)<-id2phen4(colnames(newinput))
input<-data.frame(newinput,GSE41826$beta,GSE89707$beta,GSE74486$beta,GSE103659$beta,GSE66351$beta,GSE114534$beta,check.names = F)
newinput<-input[match(rownames(newphen),colnames(input))]
GBM<-list()
GBM$beta<-newinput
GBM$phen<-newphen
save(GBM,file="MCRI.GBM.full.RData")


newphen<-subset(phen,tissue=="brain")
colnames(newinput)<-id2phen4(colnames(newinput))
input<-data.frame(newinput,GSE41826$beta,GSE89707$beta,GSE74486$beta,GSE103659$beta,GSE66351$beta,GSE114534$beta,check.names = F)
newinput<-input[match(rownames(newphen),colnames(input))]
GBM<-list()
GBM$beta<-newinput
GBM$phen<-newphen
save(GBM,file="MCRI.GBM.brain.RData")

load("MCRI.GBM.RData")
load("MCRI.GBM.full.RData")

beta<-GBM$beta
phen<-GBM$phen
beta[1:5,1:5]
phen[1:5,1:3]
phen$disease=as.character(phen$disease)
phen$tissue=as.character(phen$tissue)
phen$gender=phen$gender+1
phen$disease[phen$disease=="Normal"]<-0
phen$disease[grep("GBM|LGG",phen$disease)]<-1
phen<-data.frame(data.matrix(phen[,1:3]))
beta<-RawNARemove(beta)
beta<-data.matrix(beta)
beta<-na.omit(beta)


# 30% missing ratio, 384781 probe left. 
rlt<-c()
for(i in 1:nrow(beta)){
temp<-summary(glm(phen$disease~beta[i,]+phen$age+phen$gender,family = "binomial"))
rlt<-rbind(rlt,temp$coefficients[2,])
print(i)
}
rownames(rlt)<-rownames(beta)
fdr=p.adjust(rlt[,4],method="fdr")
newrlt<-data.frame(rlt,fdr)
save(newrlt,file="rlt.RData")
write.table(newrlt,file="MCRI.GBM.GLM.adjust.age.gender.Pval.txt",sep="\t",quote=F,col.names = NA,row.names = T)


pca <- prcomp(t(beta),center=T,scale = F)
pdf("MCRI.GBM.BUR.PCA_SDEV.pdf")
plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
dev.off()
scores <- data.frame(phen, pca$x[,1:10])
pdf("MCRI.GBM.BUR.PCA_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()
pdf("MCRI.GBM.BUR.PCA_2_3.pdf")
plot(x=scores$PC2,y=scores$PC3, xlim=c(min(scores$PC2),max(scores$PC2)),ylim=c(min(scores$PC3),max(scores$PC3)),xlab="PC2",ylab="PC3",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()
pdf("MCRI.GBM.BUR.PCA_2_4.pdf")
plot(x=scores$PC2,y=scores$PC4, xlim=c(min(scores$PC2),max(scores$PC2)),ylim=c(min(scores$PC4),max(scores$PC4)),xlab="PC2",ylab="PC4",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()
pdf("MCRI.GBM.BUR.PCA_3_4.pdf")
plot(x=scores$PC3,y=scores$PC4, xlim=c(min(scores$PC3),max(scores$PC3)),ylim=c(min(scores$PC4),max(scores$PC4)),xlab="PC3",ylab="PC4",pch=16,col=as.numeric(as.factor(phen$disease)))
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()

load("~/hpc/methylation/Normal.PBMC.GEO.HM450K.Beta.RData")
BUR<-subset(normalpbmc450beta,mean<0.1 & Quantile.75.<0.1)
marker<-read.table(file="MCRI.GBM.GLM.adjust.age.gender.Pval.txt",sep="\t",head=T,row.names = 1)
best<-marker[head(order(marker$fdr,decreasing=F),n=1000),]
best<-rownames(best[na.omit(match(rownames(BUR),rownames(best))),])
newbeta<-beta[match(best,rownames(beta)),]

GBMBUR<-list()
GBMBUR$beta<-beta
GBMBUR$phen<-phen
GBMBUR$pca
save(GBMBUR,file="MRCI.GBM.BUR.RData")

pca <- prcomp(t(beta),center=T,scale = F)

RG1<-rotation[head(order(rotation[,1],decreasing=T),n=10),]
RG2<-rotation[head(order(rotation[,1],decreasing=T),n=10),]
RG<-c(rownames(RG1),rownames(RG1))
heatmapinput<-beta[match(RG,rownames(beta)),]
mydata=t(heatmapinput)
phen=phen$disease
prefix="MCRI.GBM"
library("tsne")
library("ggplot2")
colors = rainbow(length(unique(phen)))
names(colors) = unique(phen)
ecb = function(x,y){ plot(x,t='n',xlab="Coordinate 1",ylab="Coordinate 2"); text(x,labels=phen, col=colors[phen]) }
tsne_rlt = data.frame(tsne(mydata, epoch_callback = ecb, perplexity=50))
colnames(tsne_rlt)<-c("xtsne","ytsne")
chart = ggplot(data.frame(tsne_rlt), aes(xtsne,ytsne))+geom_point(size=1,alpha=1,aes(colour = phen))+ggtitle("tSNE dimensions colored by digit")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(file="Tsne-20PC1-20PC2.pdf")

################################################################################################################
################################################################################################################
################################################################################################################
rlt<-Table2Generator(newdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newrlt<-data.frame(map[match(rownames(rlt),map[,4]),c(1:5)],rlt)
write.table(newrlt,file="../TCGA_BRCA_meth450_marker.txt",col.names=T,row.names = F,quote=F,sep="\t")
BRCAsig<-subset(newrlt,Pvalue<10^-15)
dim(BRCAsig)
head(colnames(input))
cas<-grep("Tumor",phe[,2])
con<-grep("_Normal",phe[,2])
pvalue<-apply(data,1,function(x) t.test(x[cas],x[con],paired=T)$p.value)
tvalue<-apply(data,1,function(x) t.test(x[cas],x[con],paired=T)$statistic)
output<-data.frame(data,tvalue,pvalue)
output<-subset(output,pvalue<0.0005)
dim(output)
write.table(output,file="CCA.lncRNA.diff.txt",sep="\t",quote=F,col.names = NA,row.names = T)
pdf("qqplot.MCRI.GBM.BUR.wilcox.pdf")
pQQ(marker[,4], nlabs =nrow(marker), conf = 0.95) 
dev.off()

Pvalue=apply(beta,1,function(x) wilcox.test(x[seq.control], x[seq.case],na.rm=T)$p.value)



rlt<-Table2Generator(GBM)
methdata<-GBM

Table2Generator = function(methdata){
  beta<-methdata$beta  
  phen<-methdata$phen
  beta<-data.matrix(beta)
  beta<-RawNARemove(beta)
  seq.case = which(phen$disease==1)
  seq.control = which(phen$disease == 0)
  rlt<-c()
  for(i in 1:nrow(beta)){
  temp<-t.test(beta[i,]~phen$disease)
  rlt<-rbind(rlt,c(temp$estimate,temp$statistic,temp$p.value))
  }
  colnames(rlt)<-c("Normal","Case","statistic","pval")
  Pvalue=p.adjust(rlt$pval,method="fdr")
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:nrow(beta)){
    glm.fit  = glm(phen$disease~beta[i,]+phen$age+phen$gender,family = "binomial")
    OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
    Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
    CI.upper[i]=log(exp(confint(glm.fit)[2,2]),base = 10)
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
    predicted.value = predict(glm.fit)
    predicted.data  = data.frame(Type=phen$disease, predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]]
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = logistic.rocdata[seq.max,1]
    Spec[i] = logistic.rocdata[seq.max,2]
    print(i)
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}



################################################################################################################
################################################################################################################
################################################################################################################
TSNEAnalysis<-function(mydata,phen,prefix="TSNE"){
  print("Missing Value: NA is not allowed in the input matrix!")
  print("Please be sure row is individual, column is variable!")
  print("t-sne analysis: N rows (objects) x P columns (variables) [same as PCA], be sure rownames should be unique")
  library("tsne")
  library("ggplot2")
  mydata<-na.omit(mydata)
  output=paste(prefix,".tsne-1.pdf",sep="")
  pdf(output)
  colors = rainbow(length(unique(phen)))
  names(colors) = unique(phen)
  ecb = function(x,y){ plot(x,t='n',xlab="Coordinate 1",ylab="Coordinate 2"); text(x,labels=phen, col=colors[phen]) }
  tsne_rlt = data.frame(tsne(mydata, epoch_callback = ecb, perplexity=50))
  dev.off()
  colnames(tsne_rlt)<-c("xtsne","ytsne")
  output=paste(prefix,".tsne.pdf",sep="")
  pdf(output)
  chart = ggplot(data.frame(tsne_rlt), aes(xtsne,ytsne))+geom_point(size=1,alpha=1,aes(colour = factor(phen)))+ggtitle("tSNE dimensions colored by digit")
  chart
  dev.off()
  
}

################################################################################################################
################################################################################################################
################################################################################################################

load("~/hpc/methylation/Normal.PBMC.GEO.HM450K.Beta.RData")
BUR<-subset(normalpbmc450beta,mean<0.1 & Quantile.75.<0.1)
marker<-read.table(file="MCRI.GBM.GLM.adjust.age.gender.Pval.txt",sep="\t",head=T,row.names = 1)
best<-marker[head(order(marker$fdr,decreasing=F),n=5000),]
best<-rownames(best[na.omit(match(rownames(BUR),rownames(best))),])
newbeta<-beta[match(best,rownames(beta)),]
colnames(newbeta)<-phen$disease
save(newbeta,file="newbeta.heatmap.RData")

load("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/GBM/GEO/newbeta.heatmap.RData")
HeatMap(na.omit(newbeta))


HeatMap<-function(data){
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  Hclust=function(x){hclust(x,method="complete")}
  Distfun=function(x){as.dist((1 - cor(t(x),method = "pearson")))}
  ColSideColors=sidecol(colnames(data))
  heatmap.2(data,trace="none", hclust=Hclust,
            distfun=Distfun, 
            cexRow = 1, cexCol = 1,scale="row",
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=T,Rowv = TRUE,
            keysize=0.9, margins = c(5, 10)
  )
}

################################################################################################################
################################################################################################################
################################################################################################################
load("MCRI.GBM.full.RData")
beta<-GBM$beta
phen<-GBM$phen
beta[1:5,1:5]
phen[1:5,1:3]
phen$disease=as.character(phen$disease)
phen$tissue=as.character(phen$tissue)

pca <- prcomp(t(beta),center=T,scale = F)
scores <- data.frame(phen, pca$x[,1:10])
# PCA to full dataset with blood and whole probes
pdf("MCRI.GBM.PCA.full.dataset_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$dataset))+1)
phen$col=as.numeric(as.factor(phen$dataset))+1
legend("topright",legend=c("GSE103659","GSE114534","GSE41826","GSE66351","GSE74486","GSE89707","TCGA-Brain"),pch=16,col=2:8,bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.PCA.full.tissue_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$tissue))+1)
phen$col=as.numeric(as.factor(phen$tissue))+1
legend("topright",legend=c("blood","brain"),pch=16,col=2:3,bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.PCA.full.gender_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$gender))+1)
phen$col=as.numeric(as.factor(phen$gender))
legend("topright",legend=c("male","female"),pch=16,col=c(2,3),bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.PCA.full.disease_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:4,bty="n",cex=1)
dev.off()

# PCA to full dataset with blood and BUR probes
load("~/hpc/methylation/Normal.PBMC.GEO.HM450K.Beta.RData")
BUR<-subset(normalpbmc450beta,mean<0.1 & Quantile.75.<0.1)
beta<-beta[match(rownames(BUR),rownames(beta)),]
beta<-na.omit(beta)
pca <- prcomp(t(beta),center=T,scale = F)
scores <- data.frame(phen, pca$x[,1:10])

pdf("MCRI.GBM.BUR.PCA.full.dataset_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$dataset))+1)
phen$col=as.numeric(as.factor(phen$dataset))+1
legend("topright",legend=c("GSE103659","GSE114534","GSE41826","GSE66351","GSE74486","GSE89707","TCGA-Brain"),pch=16,col=2:8,bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.BUR.PCA.full.tissue_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$tissue))+1)
phen$col=as.numeric(as.factor(phen$tissue))+1
legend("topright",legend=c("blood","brain"),pch=16,col=2:3,bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.BUR.PCA.full.gender_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$gender))+1)
phen$col=as.numeric(as.factor(phen$gender))
legend("topright",legend=c("male","female"),pch=16,col=c(2,3),bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.BUR.PCA.full.disease_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:4,bty="n",cex=1)
dev.off()
################################################################################################################
################################################################################################################
################################################################################################################
load("MCRI.GBM.brain.RData")
beta<-GBM$beta
phen<-GBM$phen
beta[1:5,1:5]
phen[1:5,1:3]
phen$disease=as.character(phen$disease)
phen$tissue=as.character(phen$tissue)
phen$disease[phen$disease=="Normal"]<-0
phen$disease[grep("GBM|LGG",phen$disease)]<-1
head(phen)
load("~/hpc/methylation/Normal.PBMC.GEO.HM450K.Beta.RData")
BUR<-subset(normalpbmc450beta,mean<0.1 & Quantile.75.<0.1)
beta<-beta[match(rownames(BUR),rownames(beta)),]
beta<-data.matrix(beta)

beta<-RawNARemove(beta)
P<-c()
for(i in 1:nrow(beta)){
  p<-t.test(beta[i,]~phen$disease)$p.value
  P<-c(P,p)
  print(i)
}

Mean<-t(apply(beta,1,function(x) tapply(x,phen$disease,function(x) mean(x,na.rm=T))))
DIFF<-Mean[,2]-Mean[,1]
rlt<-data.frame(Mean,DIFF,P)

J<-c()
NUM<-c()
for(j in seq(0.2,0.4,by=0.05)){
NUM<-c(NUM,nrow(subset(rlt,DIFF>j)))
J<-c(J,j)
}
names(NUM)<-J
pdf("markerNumber.pdf")
barplot(NUM,col="red")
dev.off()

rlt<-data.frame(Mean,DIFF,P)
rlt<-subset(rlt,DIFF>0.2)
BUR2<-read.table("~/hpc/db/hg19/GPL13534.BUR.GRCH37.merge.hg19.bed.overlap.uni.bed",sep="\t")
Tfbs<-read.table("~/hpc/db/hg19/GPL13534.wgEncodeRegTfbsClusteredV3.merge.hg19.bed.overlap.bed",sep="\t")
rlt<-rlt[rownames(rlt) %in% BUR2[,4],]
rlt<-rlt[rownames(rlt) %in% Tfbs[,4],]

for i in `ls *hg19.bed.overlap.bed`
do
sort -u $i > $i.uni.bed
done

file=list.files(pattern="*overlap.bed.uni.bed")
CPG<-c()
for(i in 1:length(file)){
cpg<-as.character(unique(read.table(file[i],sep="\t")[,4]))
CPG<-c(CPG,cpg)
}
LOCI<-read.table(file="~/hpc/db/hg19/GPL13534.RegulatoryLoci.txt",sep="\t",head=T,row.names = 1)
LOCI<-subset(LOCI,Freq==6)
rlt<-rlt[rownames(rlt)%in% rownames(LOCI),]

HeatMap(na.omit(newbeta))

newbeta<-beta[rownames(beta)%in%rownames(rlt),]

save(newbeta,file="beta.3385.RData")
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/GBM/GEO/beta.3385.RData")

colnames(newbeta)<-phen$disease
HeatMap(na.omit(newbeta))

HeatMap<-function(data){
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  Hclust=function(x){hclust(x,method="complete")}
  Distfun=function(x){as.dist((1 - cor(t(x),method = "pearson")))}
  ColSideColors=sidecol(colnames(data))
  heatmap.2(data,trace="none", hclust=Hclust,
            distfun=Distfun, 
            cexRow = 1, cexCol = 1,scale="row",
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=T,Rowv = TRUE,
            keysize=0.9, margins = c(5, 10)
  )
}



tpval<-apply(beta,1,function(x) t.test(x ~ phen$disease,na.rm=T)$p.value)
summary(glm(phen$disease~beta[1,]+ph))
Median<-t(apply(beta,1,function(x) tapply(x,phen$disease,function(x) median(x,na.rm=T))))
lapply(beta,function(x) summary(aov(x ~ phen$disease)))
lapply(beta,function(x) glm(x ~ phen$disease))

library("ReactomePA")
library("org.Hs.eg.db")
library("ggplot2")

newrlt<-data.frame(rlt,GPL[match(rownames(rlt),GPL[,4]),])
save(newrlt,file="MCRI.GBM.3385CpG.1389Gene.RData")
GPL<-read.table("~/hpc/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
GENE<-unique(na.omit(GPL[match(rownames(rlt),GPL[,4]),5]))
write.table(GPL[match(rownames(rlt),GPL[,4]),],file="MCRI.GBM.3385CpG.1389Gene.txt",sep="\t",quote=F,row.names = F,col.names = F)

load("MCRI.GBM.3385CpG.1389Gene.RData")
symbols=as.character(unique(na.omit(newrlt$V5)))
ENTREZID<-mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
x <- enrichPathway(gene=ENTREZID,pvalueCutoff=0.05, readable=T)
barplot(x, showCategory=30)
ggsave(file="MCRI.GBM.Reactome.barplot.pdf")
dotplot(x, showCategory=30)
ggsave(file="MCRI.GBM.Reactome.dotplot.pdf")
emapplot(x)
ggsave(file="MCRI.GBM.Reactome.emap.plot.pdf")
cnetplot(x, categorySize="pvalue", foldChange=geneList)
ggsave(file="MCRI.GBM.Reactome.cnet.plot.pdf")

load("continoue.RData")

LGG<-subset(phen,disease=="LGG")
GBM<-subset(phen,disease=="GBM")
Normal<-subset(phen,disease=="Normal")

hist(LGG$age)
hist(GBM$age)
hist(Normal$age)

boxplot(age~disease,data=phen,col=c(2,3,4),cex=1,pch=16,cex.axis=1.5,ylab="Age",cex.label=1.5)
summary(aov(age~disease,data=phen))

barplot(gender~disease,data=phen,col=c(2,3,4),cex=1,pch=16,cex.axis=1.5,ylab="Age",cex.label=1.5)

par(cex.axis=1.5,cex.lab=1.5)
boxplot(age~disease,data=phen,col=c(2,3,4),cex=1,pch=16,cex.axis=1.5,ylab="Age",cex.label=1.5)

par(cex.axis=1.5,cex.lab=1.5)
barplot(table(phen$gender,phen$disease),col=2:3,ylab="Sample Size",ylim=c(0,600))
legend("topleft",inset=0.05,legend = c("Male","Female"),fill=2:3,cex=1.5)


gender<-phen$gender
gender[is.na(gender)]<-"Missing"
gender[gender==0]<-"Male"
gender[gender==1]<-"Female"
pie(table(gender),col=2:4,cex=1.5)

par(cex.axis=1.25,cex.lab=1.25)
load("NUM.RData")
barplot(NUM,col=2:7,ylab="Number of DNA methylation Biomarker")


load("MCRI.GBM.MachineLearning.RData")
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/GBM/GEO/beta.3385.RData")
data<-read.table("MCRI.GBM.3358.marker.txt")
beta<-GBM$beta
phen<-GBM$phen
beta<-na.omit(beta[match(rownames(data),rownames(beta)),])
phen<-as.character(phen$disease)
data<-data.frame(y,t(beta))
model <- svm(y~ ., data = data)
pred <- predict(model, t(beta))
pred <- fitted(model)
table(pred, y)
pred <- predict(model, t(beta), decision.values = TRUE)
attr(pred, "decision.values")[1:4,]
model<-randomForest(phen ~ ., data = data, ntree = 500, mtry = 5, importance = TRUE)
load("MCRI.GBM.RandomForest.RData")
importance<-model$importance[,5]
importance<-importance[order(importance,decreasing = T)]
importance[1:50]
for(i in importance[1:50]){
  temp=data.frame(pheno=phen$disease,t(beta[match(i,rownames(beta)),]))
  model<-randomForest(pheno ~ ., data = temp, ntree = 500, mtry = 2, importance = TRUE)
  
}
load("error.RData")
par(cex.lab=1.5,cex.axis=1.5)
barplot(error,col=1:20,ylim=c(0,0.05),ylab="OOB error rate",xlab="Number of variables in the prediction model")


library(randomForest)
load("MCRI.GBM.MachineLearning.RData")
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/GBM/GEO/beta.3385.RData")
data<-read.table("MCRI.GBM.3358.marker.txt")
beta<-GBM$beta
phen<-GBM$phen
beta<-na.omit(beta[match(rownames(data),rownames(beta)),])
data<-data.frame(phen=phen$disease,t(beta))
model<-randomForest(phen ~ ., data = data, ntree = 500, mtry = 1, importance = TRUE)
write.table(model$importance,file="MCRI.GBM.RandomForest.VariableImportance.MultiClass.txt",sep="\t",quote=F,col.names = NA,row.names = T)
importance<-model$importance[,5]
importance<-importance[order(importance,decreasing = T)]
importance[1:50]

error<-c()
i=1
temp=data.frame(pheno=phen$disease,beta[match(names(importance)[1:i],rownames(beta)),])
model<-randomForest(pheno ~ ., data = temp, ntree = 500, mtry = 5, importance = TRUE)
error<-c(error,mean(model$err.rate[,1]))
for(i in 2:20){
  temp=data.frame(pheno=phen$disease,t(beta[match(names(importance)[1:i],rownames(beta)),]))
  model<-randomForest(pheno ~ ., data = temp, ntree = 500, mtry = 2, importance = TRUE)
  error<-c(error,mean(model$err.rate[,1]))
  print(i)
}

load("error.RData")
par(cex.lab=1.5,cex.axis=1.5)
barplot(error,col=1:20,ylim=c(0,0.05),ylab="OOB error rate",xlab="Number of variables in the prediction model")

load("error.3.RData")
par(cex.lab=1.5,cex.axis=1.5)
barplot(error,col=1:20,ylim=c(0,0.4),ylab="OOB error rate",xlab="Number of variables in the prediction model")



load("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/MCRI.GBM.OS.HR.1388.RData")
head(newrlt)
library("Haplin")
pdf("qqplot.MCRI.GBM.BUR.wilcox.pdf")
pQQ(newrlt[,5], nlabs =nrow(newrlt), conf = 0.95) 
dev.off()

xx<-newrlt[order(newrlt[,5],decreasing = F),]
head(xx)


load("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/MCRI.GBM.OS.HR.395788.RData")
head(rlt)
library("Haplin")
pdf("qqplot.MCRI.GBM.BUR.wilcox.pdf")
pQQ(rlt[,5], nlabs =nrow(rlt), conf = 0.95) 
dev.off()
xx<-rlt[order(rlt[,5],decreasing = F),]
head(xx)

plot(KM0, main=expression(paste("Kaplan-Meier-estimate ", hat(Lambda)(t))),xlab="t",, fun="cumhaz", lwd=2)

library(survival)
i<-match("cg04866162",rownames(newdata))
dat<-data.frame(Rna=newdata[i,],phen)
dat$Rna[dat$Rna<=0.3]<-0
dat$Rna[dat$Rna>0.3]<-1

KM0 <- survfit(Surv(month,censored)~1,  type="kaplan-meier", conf.type="log", data=dat)
KM <- survfit(Surv(month,censored)~Rna, type="kaplan-meier", conf.type="log", data=dat)

HR<-list()
HR$dat=dat
HR$KM0<-KM0
HR$KM<-KM
save(HR,file="/home/guosa/hpc/methylation/GBM/GEO/cg04866162.HR.OS.RData")
load("cg04866162.HR.OS.RData")

dat<-HR$dat
KM<-HR$KM
KM0<-HR$KM0

plot(HR$KM0,xlab="t", ylab="Survival", lwd=2)
pdf("/home/guosa/hpc/methylation/GBM/HSPB1.HR.OS.pdf")
plot(KM0, xlab="Month", ylab="Survival", lwd=2, col=1)
par(cex.lab=1.5,cex.axis=1.5)
plot(KM, xlab="Month", ylab="Survival", lwd=4, col=2:3,cex=1.5)
legend(x="topright", col=2:3, lwd=3, legend=c("Hyper","Hypo"),bty="n",cex=1.5)
dev.off()



load("MCRI.GBM.MachineLearning.RData")
load("beta.3385.RData")
data<-read.table("MCRI.GBM.3358.marker.txt")
beta<-GBM$beta
phen<-GBM$phen
beta<-na.omit(beta[match(rownames(data),rownames(beta)),])
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
rlt<-data.frame(data,map[match(rownames(data),map[,4]),])
GBMOS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/gbm/prognostic_glioma_Favourable_68_Gene.txt",sep="\t")

newrlt<-na.omit(rlt[rlt$V5 %in% GBMOS$V1,])
subset(newrlt,coef>0 & newrlt[,5]<0.005)

GBMRLT<-rlt[order(rlt[,5],decreasing = F),]


newrlt<-na.omit(rlt[rownames(rlt) %in% rownames(data),])
pQQ(newrlt[,5], nlabs =nrow(newrlt), conf = 0.95) 


head(newrlt)
load("OS.HR.LGG.RData")


GBM<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/TCGA_HM450_GBM_Suvival_HR.txt",head=T,sep="\t",row.names = 1)
LGG<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/TCGA_HM450_LGG_Suvival_HR.txt",head=T,sep="\t",row.names = 1)
data<-read.table("MCRI.GBM.3358.marker.txt")
GBM<-na.omit(GBM)
LGG<-na.omit(LGG)
head(GBM)
head(LGG)
pQQ(GBM[sample(1:153431,3000),5], nlabs =3000, conf = 0.95) 
