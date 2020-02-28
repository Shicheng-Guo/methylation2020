
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/methdata.pancancer.RData")
methdata[1:5,1:5]
methdata<-methdata[,grep("LGG",colnames(methdata))]
phen4<-id2phen4(colnames(methdata))
phen3<-id2phen3(colnames(methdata))
bin<-id2bin(colnames(methdata))
pid<-id2pid(colnames(methdata))
phen<-data.frame(phen4=phen4,phen3=phen3,pid=pid,bin=bin)
exclude<-which(c(phen$bin !=1))
phen<-phen[-exclude,]
input<-methdata[,-exclude]
Seq<-paste(phen$pid,phen$bin,sep="-")
head(phen)
input[1:5,1:5]

LGG<-grep("LGG",colnames(input))
newinput<-input[,LGG]
newphen<-phen[LGG,]
newinput[1:5,1:5]
library("survival")
library("survminer")
OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
data<-newinput[,which(id2bin(colnames(newinput))==1)]
newdata<-data[,id2phen3(colnames(data)) %in% OS$submitter_id]
colnames(newdata)<-id2phen3(colnames(newdata))
newdata<-RawNARemove(newdata)
phen<-OS[match(colnames(newdata),OS$submitter_id),]
head(phen)
phen$censored<-as.numeric(! phen$censored)
phen$month=phen$time/30
head(phen)

data<-read.table("MCRI.GBM.3024.binary.RandomForest.ImportanceVariable.txt")
newdata<-newdata[rownames(newdata)%in%rownames(data),]
HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=0.3]<-0
  dat$Rna[dat$Rna>0.3]<-1
  hr<-summary(coxph(Surv(month,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
rownames(HR)<-rownames(newdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
rlt<-data.frame(HR,map[match(rownames(HR),map[,4]),])
write.table(rlt,file="~/hpc/methylation/TCGA_HM450_LGG_Suvival_HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)

newHR<-subset(HR,HR[,5]<10^-10)
library(survival)
i<-match("cg00026222",rownames(newdata))
dat<-data.frame(Rna=newdata[i,],phen)
dat$Rna[dat$Rna<=0.3]<-0
dat$Rna[dat$Rna>0.3]<-1
KM <- survfit(Surv(month,censored)~Rna, type="kaplan-meier", conf.type="log", data=dat)
HR<-list()
HR$dat=dat
HR$KM<-KM
save(HR,file="cg00026222.HR.OS.RData")

load("cg00026222.HR.OS.RData")

dat<-HR$dat
KM<-HR$KM

pdf("/home/guosa/hpc/methylation/GBM/HSPB1.HR.OS.pdf")
plot(KM0, xlab="Month", ylab="Survival", lwd=2, col=1)
par(cex.lab=1.5,cex.axis=1.5)
plot(KM, xlab="Month", ylab="Survival", lwd=4, col=2:3,cex=1.5)
legend(x="topright", col=2:3, lwd=3, legend=c("Hyper","Hypo"),bty="n",cex=1.5)
dev.off()



rm(list = ls())
library("randomForest")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/GBM/GEO")
load("MCRI.GBM.3385CpG.1389Gene.RData")
load("MCRI.GBM.MachineLearning.RData")
beta<-GBM$beta
phen<-GBM$phen
beta<-na.omit(beta[match(rownames(data),rownames(beta)),])
phen<-as.character(phen$disease)
phen[phen=="LGG"]<-"GBM"
data<-data.frame(phen,t(beta))
model<-randomForest(phen ~ ., data = data, ntree = 500, mtry = 5, importance = TRUE)
head(model$importance)
importance<-model$importance[order(model$importance[,4],decreasing = T),]
head(importance)
write.table(importance[1:250,],file="MCRI.GBM.3024.binary.RandomForest.ImportanceVariable.txt",col.names = NA,row.names = T,sep="\t",quote=F)

beta<-beta[match(rownames(importance)[1:250],rownames(beta)),]

cpg<-"cg12974637"
x<-beta[match(cpg,rownames(beta)),]
par(cex.lab=1.5,cex.axis=1.5)
boxplot(x~phen$disease,col=2:3,main="CCDC68(cg12974637)",ylab="DNA methylation level (beta)")

phen$disease[phen$disease=="LGG"]<-"GBM"
G<-which(phen$disease=="GBM")
N<-which(phen$disease=="Normal")
N=30
library("arm")

POWER<-c()
for(z in 1:30){
power<-c()
for(i in seq(5,N,by=1)){
  pow<-c()
  for(j in 1:100){
    j1<-sample(G,i,replace = T)
    j2<-sample(N,i,replace = T)
    dd<-beta[,c(j1,j2)]
    pp<-phen[c(j1,j2),]
    glm.rlt <- bayesglm(as.factor(pp$disease)~dd[50,]+pp$gender+pp$age,family="binomial")
    fit=summary(glm.rlt)
    pow<-c(pow,fit$coefficients[2,4])
    }
    power<-rbind(power,pow)
    }
    POWER<-rbind(POWER,unlist(apply(power,1,function(x) sum(x<0.05))))
    print(z)
}

minsam<-unlist(lapply(apply(POWER,1,function(x) which(x>80)),function(x) x[1]))


plot(sort(minsam),col=1:30,cex=2,pch=16,ylab="Sample Size", xlab="Biomarker Index",)

names(powers)<-seq(5,N,by=1)
plot(powers/100,ylab="power",xlab="Sample Size",col=1:30,pch=16,cex=2)
abline(h=0.8,lwd=2,lty=2,col="blue")
abline(v=12,lwd=2,lty=2,col="blue")

rownames(FIT)<-seq(5,N,by=1)
par(cex.lab=1.5,cex.axis=1.5)
library("Haplin")
pQQ(FIT[,4],nlabs =nrow(FIT), conf = 0.95) 
plot(y=-log(FIT[,4],10),x=seq(5,N,by=1),pch=16,ylab="-log(Pval,10)",xlab="Sample Size",col=1:length(FIT[,4]),cex=2)
abline(h=-log(0.05,10),lty=2,lwd=2,col="blue")
abline(h=-log(0.001,10),lty=2,lwd=2,col="blue")




for(i in importance[1:50]){
  temp=data.frame(pheno=phen$disease,t(beta[match(i,rownames(beta)),]))
  model<-randomForest(pheno ~ ., data = temp, ntree = 500, mtry = 2, importance = TRUE)
  
}
load("error.RData")
par(cex.lab=1.5,cex.axis=1.5)
barplot(error,col=1:20,ylim=c(0,0.05),ylab="OOB error rate",xlab="Number of variables in the prediction model")
