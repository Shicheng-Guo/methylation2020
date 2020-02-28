source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_mh450/meth450Pancancer.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/meth450Pancancer.R")
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
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
NAN<-apply(input,1,function(x) any(is.na(x)))
input<-input[-which(unlist(NAN)),]
save(input,file="methdata.pancancer.nomissing.RData")


newinput<-input[,which(phen$bin==11)]
newphen<-phen[which(phen$bin==11),]

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.nomissing.RData")

random<-sample(1:ncol(input),200)
newbeta<-input[,random]
newphen3<-phen3[random]
save<-list()
save$newbeta<-newbeta
save$newphen3<-newphen3
save(save,file="randomSave.RData")
load("randomSave.RData")

beta<-save$newbeta

clinical<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/clinical.tsv",head=T,sep="\t")
clinical$age=2019-as.numeric(as.character(clinical$year_of_birth))
clinical$gender=clinical$gender
clinical$tissue=clinical$tissue_or_organ_of_origin
phen4<-id2phen4(colnames(beta))
phen3<-id2phen3(colnames(beta))
bid<-id2bin(colnames(beta))
pid<-id2pid(colnames(beta))
age<-clinical[match(phen3,clinical$submitter_id),]$age
gender<-clinical[match(phen3,clinical$submitter_id),]$gender
disease=paste(pid,bid,sep="_")
phen<-data.frame(age,disease,gender)
rownames(phen)<-phen4

pca <- prcomp(t(newinput),center=F,scale = F)
pdf("MCRI.GBM.BUR.PCA_SDEV.pdf")
plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
dev.off()

scores <- data.frame(newphen, pca$x[,1:10])
pdf("MCRI.TCGA.Normal.PCA_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=names(table(phen$disease)),pch=16,col=1:length(names(table(phen$disease))),bty="n",cex=0.5)
dev.off()

scores <- data.frame(newphen, pca$x[,1:10])
pdf("../MCRI.TCGA.Normal.PCA_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$pid))+1)
phen$col=as.numeric(as.factor(phen$pid))+1
legend("topright",legend=names(table(phen$pid)),pch=16,col=1:length(names(table(phen$pid))),bty="n",cex=0.5)
dev.off()

d <- dist(t(newinput)) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="Metric MDS", type="n")
text(x, y, labels = row.names(mydata), cex=.7)






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



