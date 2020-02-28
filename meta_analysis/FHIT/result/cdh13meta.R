# pipeline meta analysis for CHD13 methylation diagnostic role in NSCLC. 
setwd("C:/Users/shicheng/Dropbox/Project/methylation/meta/FHIT")
install.packages("meta")
install.packages("metafor")
install.packages("mada")

library("meta")
library("metafor")
library("mada")

apc<-read.table("FHIT.txt",head=T,sep="\t")
head(apc)

# 1, meta analysis and forest plot
meta1 <- metabin(event.e, n.e, event.c, n.c,data=apc,studlab=author,sm="OR", method="I",method.tau="DL")
meta1
pdf("Figure1.forest.pdf")
forest.meta(meta1,fontsize=6,leftlabs=c(NA, "Methylated", NA, "Methylated", NA),lab.e="NSCLC",lab.c="Normal",rightcols=c("effect","ci","w.random"),comb.random=TRUE,comb.fixed=TRUE,cex=0.4)
dev.off()

# 2, cumulative analysis by the order of publish year
pdf("Figure2.cum_year.pdf")
cum<-metacum(meta1,sortvar=apc$year,pooled="random")
forest.meta(cum,fontsize=6)
dev.off()

# 3, trim fill plot
tf1<-trimfill(meta1,comb.fixed=T)
tf1
pdf("Figure3.trimfill.pdf")
forest.meta(tf1,fontsize=6,comb.fixed=TRUE)  #
dev.off()

# 4, funnel plot
pdf("Figure4.funnel.pdf")
funnel.meta(meta1,log="x")
dev.off()

# 5, labbe plot
pdf("Figure5.labbe.pdf")
labbe(meta1,studlab=apc$author,cex.studlab=0.4,xlim=c(-0.2,1))
dev.off()

# 6, Egger test and Harbord test
metabias(meta1,method.bias="linreg") # Egger test for non odds ratio effect 
metabias(meta1,method.bias="score") # Harbord test for odds ratio effect

# 7,leave-one meta analysis
pdf("Figure6.leave-one-meta.pdf")
forest(metainf(meta1, pooled="random"), comb.random=TRUE,fond=4)
# forest(metainf(meta1, pooled="fixed"), comb.random=TRUE,fond=4) # 
dev.off()

# 8, combined sensitivity and specificty and SROC
TP<-apc$event.e
FN<-apc$n.e-apc$event.e
FP<-apc$event.c
TN<-apc$n.c-apc$event.c
data<-data.frame(TP,FN,FP,TN)
fit.roc<-reitsma(data)
summary(fit.roc)    # extract combined Sen and Spe from Here

tissue<-data[apc$sampletype=="Tissue",]
serum<-data[apc$sampletype=="Non-tissue",]
fit.tissue<-reitsma(tissue)
fit.serum<-reitsma(serum)
summary(fit.tissue)
summary(fit.serum)

pdf("SROC2.pdf")
plot(fit.roc, xlim = c(0,0.7), ylim = c(0,0.75),main = "Diagnostic SROC(bivariate model) for FHIT in NSCLC",cex.main=1,cex.lab=1,cex.axis=1,lwd=4,lty=1)
lines(sroc(fit.roc), lty = 2,lwd=4)
ROCellipse(fit.roc, lty = 2, pch = 2, d = TRUE)
points(fpr(data), sens(data), cex = 1,pch=2,lwd=4)
legend("bottomright", c("Observed data in NSCLC","SROC for NSCLC", "95%CI for NSCLC"), pch = c(2,NA,NA), lty = c(NA,2,1),cex=1.1,lwd=4)
dev.off()

# stepwise meta regression
subgroup<-c("sampletype","age","stage1","stage2","genderatio","methods","aim","Multipletarget","controldesign","Ad2Sc","primer")
colnames(apc)
collist<-na.omit(match(subgroup,colnames(apc)))
collist
metaregrlt<-list()

rlt1<-matrix(NA,length(collist),7)
for (i in 1:length(collist)){
  apcomit<-apc
  metaomit<-metabin(event.e, n.e, event.c, n.c,data=apcomit,sm="OR", method="I",method.tau="DL") #complab=c("NSCLC","NORMAL"),
  metaregrlt<-metareg(metaomit,~apcomit[,collist[i]]) 
  
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
rownames(rlt1)<-colnames(apc)[collist]
colnames(rlt1)<-c("Coefficient","ci.lb","ci.ub","pvalue","tau2","QE","QE.P-value")
rlt1
write.table(rlt1,file="meta-regression.txt",sep="\t",row.names=T,col.names=NA,quote=F)

# conditional meta regression
metaregrlt<-metareg(meta1,~apc$sampletype+apc$stage2+apc$aim)
metaregrlt

# stepwise stratified analysis
subgroup<-c("sampletype","age","stage1","stage2","genderatio","methods","aim","Multipletarget","controldesign","Ad2Sc","primer")
collist<-na.omit(match(subgroup,colnames(apc)))
meta1 <- metabin(event.e, n.e, event.c, n.c,data=apc,studlab=author2,sm="OR", method="I",complab=c("NSCLC","NORMAL"),method.tau="DL",byvar=apc[,collist[1]])
collist

i=3
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F)
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(4, "cm"))
dev.off()


i<-4
apc[,4]
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(subdata[,i],na.rm=T)
median
less<-which(subdata[,i]<59)
less
more<-which(subdata[,i]>=59)
more
subdata[less,i]<-"age<59"
subdata[more,i]<-"age>=59"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(5, "cm"))
dev.off()



i<-5
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=0.50)
more<-which(subdata[,i]>0.50)
subdata[less,i]<-"stage1<50.45%"
subdata[more,i]<-"stage1>=50.45%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(5, "cm"))
dev.off()

i<-6
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=0.60)
more<-which(subdata[,i]>0.60)
less
more
subdata[less,i]<-"stage <60%%"
subdata[more,i]<-"stage(I+II)>=60%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F)
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(4, "cm"))
dev.off()


i<-7
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=median)
more<-which(subdata[,i]>median)
subdata[less,i]<-"Male vs Female <71.31%"
subdata[more,i]<-"Male vs Female >=71.3%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F)
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(6, "cm"))
dev.off()

head(apc)


for (i in 12:15){
  valid<-! is.na(apc[,i])
  subdata<-apc[! is.na(apc[,i]), ]
  meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F)
  meta2
  filename<-paste(colnames(subdata)[i],"_stratify.meta.txt",sep="")
  write.table(meta2,file=filename)
  filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
  pdf(filename)
  forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(4, "cm"))
  dev.off()
}





i<-15
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
subdata[,i]<-as.array(as.character(subdata[,i]))
subdata[,i]
less<-which(subdata[,i]=="hom")
more<-which(subdata[,i]=="heter")
subdata[less,i]<-"Autogenous"
subdata[array(more),i]<-"Heterogeneous"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(4, "cm"))
dev.off()



i<-17
valid<-! is.na(apc[,i])
subdata<-apc[! is.na(apc[,i]), ]
median<-median(apc[,i],na.rm=T)
median
less<-which(subdata[,i]<=median)
more<-which(subdata[,i]>median)
less
more
subdata[less,i]<-"Ad2Sc < 79.6%"
subdata[more,i]<-"Ad2Sc >= 79.6%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=F,label.e="NSCLC",label.c="Normal")
meta2
summary(meta2)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=12,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=F,squaresize=1.6,plotwidth=grid::unit(4, "cm"))
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



