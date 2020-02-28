setwd("/home/sguo/score");
source("fclust.R")
library("skmeans")
library("plotrix")
library("knitr")
library("fpc")

# shared clinical data by seq-cluster and array-cluster
load("OvClinic.RData")
ClinicNumericFact<-matrix(NA,412,8)
for (i in 2:9){
  ClinicNumericFact[,i-1]<-as.numeric(clinc[,i])
}
colnames(ClinicNumericFact)<-tolower(colnames(clinc)[2:9])

# seq-cluster
data<-load("ScoreMatrix.RData")
SeqSkmeansRlt<-skmean(ScoreMatrix,10)  # 10 here means 2:10 not 1:10 in my kmean function;
save(SeqSkmeansRlt,file="SeqSkmeanRlt.RData")
#write.table(result,file="fhclust.txt",sep="\t",row.names=T,col.names=NA)

# array-cluster
load("ExpressDataOrdered.RData")
ArrSkmeansRlt<-skmean(Expre,10)  # 10 here means 2:10 not 1:10 in my kmean function;
save(ArrSkmeansRlt,file="ArrSkmeanRlt.RData")

# compare with clinic clust
rlt1<-ClustValidtion(ArrSkmeansRlt,ClinicNumericFact)
rlt2<-ClustValidtion(SeqSkmeansRlt,ClinicNumericFact)

# plot the result
pdf("array.cluster.pdf")
par(mfrow=c(2,2))
color2D.matplot(rlt1$vi,show.values=2,show.legend=T,ylab="Cluster By FPCA",xlab="Clinical Cluster",main="Array VI 2D matrix plot")
color2D.matplot(rlt2$vi,show.values=2,show.legend=T,ylab="Cluster By FPCA",xlab="Clinical Cluster",main="RNA-Seq VI 2D matrix plot")
color2D.matplot(rlt1$correctRand,show.values=2,show.legend=T,ylab="Cluster By FPCA",xlab="Clinical Cluster",main="Array CorrectRand 2D matrix plot")
color2D.matplot(rlt2$correctRand,show.values=2,show.legend=T,ylab="Cluster By FPCA",xlab="Clinical Cluster",main="RNA-Seq CorrectRand 2D matrix plot")
dev.off()


# compare seq result and array result(cluster validation)
rlt2<-rlt[c(which(rlt$clinc=="2"),which(rlt$clinc== "3")),]
xtable(~clinc+clust,data=rlt2)
xtabs(~clinc+clust,data=rlt2)

mydata<-matrix(rnorm(10000000,1,1),100,100000)
sparcl(mydata,method="complete",cut=4)


# sparse cluster analysis
data<-ScoreMatrix
pheno=clinc[,9]
numberFactor<-CharToNumberFactor(pheno)
data.frame(pheno,numberFactor)
dat<-ColNARemove(dat)
dat<-RawNARemove(dat)

dat<-data[c(which(pheno=="Sensitive"), which(pheno=="Resistant")),]
phen<-pheno[c(which(pheno=="Sensitive"), which(pheno=="Resistant"))]
hightaucolumn<-pvalueSparse(dat,phen,"Sensitive","Resistant",pvalue=0.005)
rlt<-sparcl(as.matrix(hightaucolumn),method="complete",cut=2)




sparsedata<-sdsparse(dat,0.75)
hightaucolumn<-sparcl(as.matrix(sparsedata),method="complete",cut=2)




# Subfunction
ClustValidtion<-function(ClusteResult,ClinicNumericFact){
rlt<-list()
data<-cbind(ClinicNumericFact,ClusteResult)
data2<-RawNARemove(data)
ClinicNumericFact<-data2[,1:dim(ClinicNumericFact)[2]]
ClusteResult<-data2[,(dim(ClinicNumericFact)[2]+1):(dim(ClinicNumericFact)[2]+dim(ClusteResult)[2])]
CorrectRand<-Vi<-matrix(NA,dim(ClusteResult)[2],dim(ClinicNumericFact)[2])
for (i in 1:dim(ClusteResult)[2]){
  for (j in 1:dim(ClinicNumericFact)[2]){
    result<-cluster.stats(NULL,ClinicNumericFact[,j],ClusteResult[,i],compareonly=T)
    CorrectRand[i,j]<-result$corrected.rand
    Vi[i,j]<-result$vi
  }
}
   rlt$correctRand<-CorrectRand
   rlt$vi<-Vi
   rlt
}

skmean<-function(data, classnum){
  # suitable to big data especially for col > row matrix, defort to cluste matrix to maximun 10 clusters
  cluster<-matrix(NA, dim(data)[1], classnum-1)
  for (i in 2:classnum){
    ask<-skmeans(data,i)
    file<-paste(data,".sk",i,".RData",sep="")
    cluster[,i-1]<-ask$cluster
    save(ask, file=file)
    print(i)
  }
  cluster
}

RawNARemove<-function(data){
  NaRaw<-which(apply(data,1,function(x) any(is.na(x)))==T)
  if(length(NaRaw)>0){
    data1<-data[-NaRaw,]
  }else{
    data1<-data;
  }
  data1
}
ColNARemove<-function(data){
  NaCol<-which(apply(data,2,function(x) any(is.na(x)))==T)
  if(length(NaCol)>0){
    data1<-data[,-NaCol]
  }else{
    data1<-data;
  }
  data1
}

sparcl<-function(mydata,method="complete",cut=4){
library("sparcl")
set.seed(1)
# Do tuning parameter selection for sparse hierarchical clustering
perm.out <- HierarchicalSparseCluster.permute(mydata, wbounds=c(1.5,2:9),nperms=5)
# Perform sparse hierarchical clustering
sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw,method="complete")
m<-mydata[,which(sparsehc$ws !=0)]
#drawHeatmap2(m)
d <- dist(m, method = "euclidean") #distance matrix
fit <- hclust(d, method="ward")
#plot(fit) # display dendogram
clusters <- cutree(fit, k=cut) #cut tree into 5 clusters
#draw dendogram with red borders around the 5 cluster
#rect.hclust(fit, k=cut, border="red") 
#retrun output
rlt<-list()
rlt$featurenumber<-sum(sparsehc$ws !=0)
rlt$group<-clusters;
rlt$sparsedata<-m
rlt
}


sdsparse<-function(dat,quantile){
  sd<-apply(dat,2,function(x) sqrt(var(x))/mean(x))
  dat<-dat[,which(sd>quantile(sd,quantile))]
  dat
}

tauSparse<-function(dat,pheno,tau=0.8){  
  library(gplots)
  data<-as.list(data.frame(dat))
  aggre<-aggregate(data,by=list(pheno),FUN=median)
  tau2<-as.numeric();
  for (i in 2:(ncol(aggre)-1)){
    tau1<-(nrow(aggre)-sum(aggre[,i])/max(aggre[,i]))/(nrow(aggre)-1)	
    tau2<-c(tau2,tau1)
  }
  
  color.map <- function(cancertype){ 
    cancertype<-as.numeric(cancertype)
    cancertype[is.na(cancertype)]<-0
    cancertype
  }
  
  patientcolors <- as.character(unlist(lapply(pheno, color.map)))
  data<-as.data.frame(data)
  pdf("heatmap.pdf",width=17,height=17)  #increase it when fig is big
  heat<-heatmap.2(as.matrix(data[,which(tau2>tau)]),col=topo.colors(75), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors)
  dev.off()
  hightaucolumn<-data[,which(tau2>tau)]
  hightaucolumn
}


CharToNumberFactor<- function(Variable){ 
  Variable<-as.numeric(Variable)
  Variable[is.na(Variable)]<-0
  Variable
}

  
pvalueSparse<-function(dat,pheno,phenocode1,phenocode2,pvalue=0.005){  
  library(gplots)
  p<-apply(dat,2, function(x) wilcox.test(x[which(pheno==phenocode1)],x[which(pheno==phenocode2)])$p.value)  
  color.map <- function(cancertype){ 
    cancertype<-as.numeric(cancertype)
    cancertype[is.na(cancertype)]<-0
    cancertype
  }
  patientcolors <- as.character(unlist(lapply(phen, color.map)))
  data<-as.data.frame(dat)
  pdf("heatmap.v=0.00005.pdf",width=17,height=17)  #increase it when fig is big
  heat<-heatmap.2(as.matrix(data[,which(p<0.005)]),col=topo.colors(75), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors,main="p-value=0.005")
  dev.off()
  hightaucolumn<-data[,which(p<pvalue)]
  hightaucolumn
}


