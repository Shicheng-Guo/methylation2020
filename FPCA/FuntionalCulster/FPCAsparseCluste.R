setwd("/home/sguo/score");
source("fclust.R")
library("skmeans")
library("plotrix")
library("knitr")
library("fpc")

# Shared clinical data by seq-cluster and array-cluster
load("OvClinic.RData")
ClinicNumericFact<-matrix(NA,412,8)
for (i in 2:9){
  ClinicNumericFact[,i-1]<-as.numeric(clinc[,i])
}
colnames(ClinicNumericFact)<-tolower(colnames(clinc)[2:9])
pheno=clinc[,9]

# seq-cluster
data<-load("ScoreMatrix.RData")
SeqSkmeansRlt<-skmean(ScoreMatrix,10)  # 10 here means 2:10 not 1:10 in my kmean function;
save(SeqSkmeansRlt,file="SeqSkmeanRlt.RData")
#write.table(result,file="fhclust.txt",sep="\t",row.names=T,col.names=NA)

# array-cluster
load("ExpressDataOrdered.RData")
# ArrSkmeansRlt<-skmean(Expre,10)  # 10 here means 2:10 not 1:10 in my kmean function;
# save(ArrSkmeansRlt,file="ArrSkmeanRlt.RData")

dat<-ExpressData2[c(which(pheno=="Sensitive"), which(pheno=="Resistant")),]
phen<-pheno[c(which(pheno=="Sensitive"), which(pheno=="Resistant"))]

rlt_1<-sparHierarchiCl(as.matrix(dat),method="complete",cut=2)

save(rlt_1,file="array_sdsparse_sparHierarchiCl.RData")
rlt_2<-sparKmeanCl(as.matrix(dat),k=2)
save(rlt_2,file="array_sdsparse_sparKmeanCl.RData")


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

# sparse cluster analysis
load("ScoreMatrix.RData")
load("OvClinic.RData")
data<-ScoreMatrix
pheno=clinc[,9]
numberFactor<-CharToNumberFactor(pheno)
data.frame(pheno,numberFactor)
data<-ColNARemove(data)
data<-RawNARemove(data)

#just for drug-response analysis
dat<-data[c(which(pheno=="Sensitive"), which(pheno=="Resistant")),]
phen<-pheno[c(which(pheno=="Sensitive"), which(pheno=="Resistant"))]
sparsedata1<-sdsparse(dat,0.25)
rlt1<-sparHierarchiCl(as.matrix(sparsedata1),method="complete",cut=i)
rlt2<-sparKmeanCl(as.matrix(sparsedata1),k=i)

save(rlt1,file="sdsparse_sparHierarchiCl.RData")
save(rlt2,file="sdsparse_sparKmeanCl.RData")


sparsedata2<-pvalueSparse(dat,phen,"Sensitive","Resistant",pvalue=0.005)
rlt3<-sparHierarchiCl(as.matrix(sparsedata2),method="complete",cut=2)
rlt4<-sparKmeanCl(as.matrix(sparsedata2),k=2)
save(rlt3,file="pvalueSparse_sparHierarchiCl.RData")
save(rlt4,file="pvalueSparse_sparKmeanCl.RData")

#comparison
cc1<-data.frame(rlt1$group,phen)
dd1<-xtabs(~cc1[,1]+cc1[,2],cc1)

cc2<-data.frame(rlt2clust,phen)
dd1<-xtabs(~cc2[,1]+cc2[,2],cc2)

cc3<-data.frame(rlt3$group,phen)
dd1<-xtabs(~cc3[,1]+cc3[,2],cc3)

cc4<-data.frame(rlt4$group,phen)
dd1<-xtabs(~cc4[,1]+cc4[,2],cc4)

cc5<-data.frame(rlt_1$group,phen)
dd1<-xtabs(~cc5[,1]+cc5[,2],cc5)

# merge cluste group 
setwd("/home/sguo/score")

load("sdsparse_sparHierarchiCl.cut=2.RData")
x2<-rlt1.2$group
load("sdsparse_sparHierarchiCl.cut=3.RData")
x3<-rlt1.3$group
load("sdsparse_sparHierarchiCl.cut=4.RData")
x4<-rlt1.4$group
load("sdsparse_sparHierarchiCl.cut=5.RData")
x5<-rlt1.5$group
load("sdsparse_sparHierarchiCl.cut=6.RData")
x6<-rlt1.6$group
load("sdsparse_sparHierarchiCl.cut=7.RData")
x7<-rlt1.7$group
load("sdsparse_sparHierarchiCl.cut=8.RData")
x8<-rlt1.8$group
load("sdsparse_sparHierarchiCl.cut=9.RData")
x9<-rlt1.9$group
load("sdsparse_sparHierarchiCl.cut=10.RData")
x10<-rlt1.10$group
clust<-cbind(x2,x3,x4,x5,x6,x7,x8,x9,x10)

ClusteResult<-clust
dim(ClinicNumericFact)
rlt<-ClustValidtion(clust,ClinicNumericFact)
write.table(rlt$correctRand,file="correctRand.txt",sep="\t",quote=F)
write.table(rlt$vi,file="vi.txt",sep="\t",quote=F)


#sumulation



#Subfunction
ClustValidtion<-function(ClusteResult,ClinicNumericFact){
  library("fpc")
  rlt<-list()
  ClusteResult<-as.matrix(ClusteResult)
  ClinicNumericFact<-as.matrix((ClinicNumericFact))
  data<-cbind(ClinicNumericFact,ClusteResult)
  data2<-RawNARemove(data)
  ClinicNumericFact<-data2[,1:dim(ClinicNumericFact)[2]]
  ClusteResult<-data2[,(dim(as.matrix(ClinicNumericFact))[2]+1):(dim(as.matrix(ClinicNumericFact))[2]+dim(ClusteResult)[2])]
  CorrectRand<-Vi<-matrix(NA,dim(ClusteResult)[2],dim(as.matrix(ClinicNumericFact))[2])
  for (i in 1:dim(ClusteResult)[2]){
    for (j in 1:dim(as.matrix(ClinicNumericFact))[2]){
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
sparHierarchiCl<-function(mydata,method="complete",cut=4){
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
  rlt$bestw<-perm.out$bestw
  rlt
}
sparKmeanCl<-function(mydata,k=2){
  library("sparcl")
  set.seed(1)
  km.perm <- KMeansSparseCluster.permute(mydata,k,wbounds=seq(3,10,len=15),nperms=5)
  # run sparse k-means
  km.out <- KMeansSparseCluster(mydata,k,wbounds=km.perm$bestw)
  rlt<-list()
  rlt$sparsedata<-mydata[,km.perm$nnonzerows]
  rlt$print<-km.out
  print(km.out,file="SparseKemansClustResult.txt")
  rlt
}
sdsparse<-function(dat,quantile){
  sd<-apply(dat,2,function(x) sqrt(var(x))/mean(x,na.rm = T))
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
    p<-apply(dat,2, function(x) wilcox.test(x[which(pheno==phenocode1)], x[which(pheno==phenocode2)])$p.value)  
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
frr<-function(rr1,rr2,p,q,lamdi,frr){
#familial reltaive risk
r1<-10   #relative risk(estimated by odds ratios) for hetrozygotes relative to common homozygotes
r2<-10   #relative risk(estimated by odds ratios) for rare homozygotes to common homozygotes
p<-0.3   #risk allelefrequency
q<-1-p
lamdi0<-8.48  #overall familial relatve risk 
frr<-(p*(p*r2+q*r1)^2+q*(p*r1+q*r2)^2)/(p^2*r2+2*p*q*r1+q^2)^2
proportion<-log(lamdi)/log(lamdi0)
frr
proportion
}
RawNARemove<-function(data,missratio=0.3){
  NaRaw<-which(apply(data,1,function(x) (is.na(x)))==T)
  if(length(NaRaw)>(1-missratio)*dim(data)[2]){
    data1<-data[-NaRaw,]
  }else{
    data1<-data;
  }
  data1
}
ColNARemove<-function(data,missratio=0.3){
  NaCol<-which(apply(data,2,function(x) (is.na(x)))==T)
  if(length(NaCol)>(1-missratio)*dim(data)[1]){
    data1<-data[,-NaCol]
  }else{
    data1<-data;
  }
  data1
}






