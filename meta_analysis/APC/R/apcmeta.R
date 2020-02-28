library("meta")
library("metafor")
library("mada")

setwd("/home/sguo/Dropbox/Project/methylation/meta/APC/R")
apc<-read.table("apcdata.txt",head=T,sep="\t")

TP<-apc$event.e
FN<-apc$n.e-apc$event.e
FP<-apc$event.c
TN<-apc$n.c-apc$event.c

data<-data.frame(TP,FN,FP,TN)
data

random<-mauni(x=data,method="DSL")
fixed<-mauni(x=data,method="MH")
forest(random,log=F)
apc$author
dim(data)

forest(random)
old.par <- par()
AuditC.d <- madad(AuditC)
print(AuditC.d, digits = 2) #round everything to 2 digits

par(old.par)

summary(random)
summary(fixed)

plot.new()
par(fig = c(0, 0.5, 0, 1), new = TRUE)
forest(AuditC.d, type = "sens", xlab = "Sensitivity")
par(fig = c(0.5, 1, 0, 1),  new = TRUE)
forest(AuditC.d, type = "spec", xlab = "Specificity")
forest(AuditC.d, type = "spec",  xlab = "Specificity",snames = LETTERS[1:14])
? madad
madad(data)
forest(madad(data),type="sens")
forest(m(data),type="spec")
ROCellipse(data, pch = "")
points(fpr(data), sens(data))
(fit.reitsma <- reitsma(data))
summary(fit.reitsma)
plot(fit.reitsma, sroclwd = 2,main = "SROC curve (bivariate model) for APC methylation data")
points(fpr(data), sens(data), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))

fit.guo<-reitsma(data)
summary(fit.guo)

fit.guo1<-reitsma(data[apc$sampletype=="serum",])
fit.guo2<-reitsma(data[apc$sampletype=="tissue",])
summary(fit.guo)
fit.guo<-reitsma(data[apc$methods=="MSP",])
summary(fit.guo)
fit.guo<-reitsma(data[apc$methods=="qMSP",])
summary(fit.guo)

pdf("SROC2.pdf")
plot(fit.guo, xlim = c(0,1), ylim = c(0,1),main = "Comparison of SROC(bivariate model) for NSCLC and Pca",cex.main=0.6,cex.lab=0.6,cex.axis=0.5)
lines(sroc(fit.yang), lty = 2)
ROCellipse(fit.yang, lty = 2, pch = 2, d = TRUE)
points(fpr(data2), sens(data2), cex = 0.6,pch=2)
points(fpr(data), sens(data), cex = 0.6, pch=1)
legend("bottomright", c("Observed data in NSCLC","Observed data in PCa","SROC for NSCLC", "SROC for PCa","95%CI for NSCLC","95%CI for PCa"), pch = c(1,2,1,2,NA,NA), lty = c(0,0,1,2,1,2),cex=0.4)
text(y=0.52,x=0.215,labels="sen=0.55,spe=0.78",cex=0.35)
text(y=0.72,x=0.155,labels="sen=0.75,spe=0.85",cex=0.35)
dev.off()
? reitsma
summary(fit.guo)




pdf("SROC2.pdf")
plot(fit.guo, xlim = c(0,1), ylim = c(0,1),main = "SROC for Non-small cell lung cancer(bivariate model)",cex.main=1,cex.lab=1,cex.axis=1)
points(fpr(data), sens(data), cex = 1, pch=1)
legend("bottomright", c("Observed data in NSCLC","SROC for NSCLC","95%CI for NSCLC"), pch = c(1,1,NA), lty = c(0,1,1),cex=1)
text(y=0.6,x=0.3,labels="sen=0.55,spe=0.78",cex=0.8,col="red")
dev.off()


? reitsma
summary(fit.guo)



              ###########################

      ######   trtional meta analysis  #########

              ###########################

library("meta")
library("metafor")
setwd("/home/gsc/Dropbox/Project/APCmeta/R")
apc<-read.table("apcdata.txt",head=T,sep="\t")
meta1 <- metabin(event.e, n.e, event.c, n.c,data=apc,studlab=author,sm="OR", method="I",method.tau="DL")
pdf("forest.eps")
forest.meta(meta1,fontsize=6,leftlabs=c(NA, "Methylated", NA, "Methylated", NA),lab.e="NSCLC",lab.c="Normal",rightcols=c("effect","ci","w.random"),comb.random=TRUE)
dev.off()



tiff("cum_year.tiff",res=300)
cum<-metacum(meta1,sortvar=apc$year,pooled="random")
forest.meta(cum,fontsize=6)
dev.off()

tf1<-trimfill(meta1)
pdf("trimfill.pdf")
forest.meta(tf1,fontsize=6)
dev.off()

pdf("funnel.pdf")
funnel.meta(meta1,log="x")
dev.off()
? funnel.meta
pdf("labbe.pdf")
labbe(meta1,studlab=apc$author,cex.studlab=0.4,xlim=c(-0.2,1))
dev.off()
metabias(meta1)
colnames(apc)

# stepwise meta regression
subgroup<-c("sampletype","age","stage1","stage2","genderatio","methods","aim","Multipletarget","controldesign","Ad2Sc","primer")
collist<-match(subgroup,colnames(apc))
metaregrlt<-list()
for (i in 1:length(collist)){
  for (j in 1:dim(apc)[1]){
    apcomit<-apc[-j,]
    metaomit<-metabin(event.e, n.e, event.c, n.c,data=apcomit,sm="OR", method="I",method.tau="DL") #complab=c("NSCLC","NORMAL"),
    metaregrlt<-metareg(~apcomit[,collist[i]], data=metaomit)  
    if(metaregrlt$pval[2]<0.05){
    combination<-paste(apc$author[j],colnames(apc)[collist[i]],metaregrlt$pval[2], sep=".")
    print(combination)
    }
  }
}

# stepwise meta regression
subgroup<-c("sampletype","age","stage1","stage2","genderatio","methods","aim","Multipletarget","controldesign","Ad2Sc","primer")
collist<-match(subgroup,colnames(apc))
metaregrlt<-list()

rlt1<-matrix(NA,11,7)
for (i in 1:length(collist)){
    apcomit<-apc
    metaomit<-metabin(event.e, n.e, event.c, n.c,data=apcomit,sm="OR", method="I",method.tau="DL") #complab=c("NSCLC","NORMAL"),
    metaregrlt<-metareg(~apcomit[,collist[i]], data=metaomit) 
    
    beta<-round(metaregrlt$b[2],2)
    ci.lb<-round(metaregrlt$ci.lb[2],2)
    ci.ub<-round(metaregrlt$ci.ub[2],2)
    pvalue<-(metaregrlt$pval[2])
    tau2<-round(metaregrlt$tau2,2)
    QE<-round(metaregrlt$QE,2)
    QEp<-(metaregrlt$QEp)
    
    combination<-paste(colnames(apc)[collist[i]],beta,"(",ci.lb,",",ci.ub,")",pvalue, sep=" ")
    print(combination)
    rlt1[i,1]<-beta
    rlt1[i,2]<-ci.lb
    rlt1[i,3]<-ci.ub
    rlt1[i,4]<-pvalue  
    rlt1[i,5]<-tau2  
    rlt1[i,6]<-QE   
    rlt1[i,7]<-QEp    
}
rownames(rlt1)<-subgroup
colnames(rlt1)<-c("Coefficient","ci.lb","ci.ub","pvalue","tau2","QE","QE.P-value")
write.table(rlt1,file="meta-regression.txt",sep="\t")

# conditional meta regression
metaregrlt<-metareg(~apc$primer+apc$age+apc$Ad2Sc+apc$controldesign, data=meta1)


# stepwise stratified analysis
subgroup<-c("sampletype","age","stage1","stage2","genderatio","methods","aim","Multipletarget","controldesign","Ad2Sc","primer")
collist<-match(subgroup,colnames(apc))
meta1 <- metabin(event.e, n.e, event.c, n.c,data=apc,studlab=author2,sm="OR", method="I",complab=c("NSCLC","NORMAL"),method.tau="DL",byvar=apc[,collist[1]])

i<-4
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<65)
more<-which(subdata[,i]>=65)
subdata[less,i]<-"age<65"
subdata[more,i]<-"age>=65"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()



i<-5
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=55)
more<-which(subdata[,i]>55)
subdata[less,i]<-"stage1<49.45%"
subdata[more,i]<-"stage1>=49.45%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()

i<-6
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=median)
more<-which(subdata[,i]>median)
subdata[less,i]<-"stage <75.33%"
subdata[more,i]<-"stage(I+II)>=75.33%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()


i<-7
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=median)
more<-which(subdata[,i]>median)
subdata[less,i]<-"Male vs Female <69.1%"
subdata[more,i]<-"Male vs Female >=69.1%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()


for (i in 12:15){
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()
}


i<-15
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
subdata[,i]<-as.array(as.character(subdata[,i]))
subdata[,i]
less<-which(subdata[,i]=="hom")
more<-which(subdata[,i]=="heter")
subdata[less,i]<-"Autogenous"
subdata[array(more),i]<-"Heterogeneous"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()



i<-17
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=1.9)
more<-which(subdata[,i]>1.9)
subdata[less,i]<-"Ad2Sc < 2"
subdata[more,i]<-"Ad2Sc >= 2"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
summary(meta2)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()




forest.meta(meta1)
pdf("metainf.eps")
forest(metainf(meta1, pooled="random"), comb.random=TRUE,fond=6)
dev.off()


i<-18
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i] =="set1")
more<-which(subdata[,i] =="set2")
subdata[array(less),i]<-"set1"
subdata[array(more),i]<-"set2"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
summary(meta2)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=unit(6, "cm"))
dev.off()




