# Quantitative assessment of the diagnostic role of CDH13 promoter methylation in non-small cell lung cancer
# Feb 15, 2015
# pipeline meta analysis for CDH13 methylation diagnostic role in NSCLC. 

source("https://bioconductor.org/biocLite.R")
biocLite("meta")
biocLite("metafor")
biocLite("mada")

library("meta")
library("metafor")
library("mada")

setwd("/home/sguo/Dropbox/Project/methylation/meta/CHD13")
setwd("C:/Users/shicheng/Dropbox/Project/methylation/meta/CHD13")
input<-read.table("cdh13.3.txt",head=T,sep="\t")
# input<-read.csv("CDH13-Independent.csv", header = T)
head(input)

# 1, meta analysis and forest plot
meta1 <- metabin(event.e, n.e, event.c, n.c,data=input,studlab=author,sm="OR", method="MH",method.tau="DL")
pdf("Figure1.forest.pdf")
forest.meta(meta1,fontsize=10,leftlabs=c(NA, "Methylated", NA, "Methylated", NA),lab.e="NSCLC",lab.c="Normal",rightcols=c("effect","ci","w.random"),comb.random=TRUE,comb.fixed=TRUE,cex=0.5)
dev.off()

# 2, cumulative analysis by the order of publish year
pdf("Figure2.cum_year.pdf")
cum<-metacum(meta1,sortvar=input$year,pooled="random")
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
labbe(meta1,studlab=input$author,cex.studlab=0.4,xlim=c(-0.2,1))
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
TP<-input$event.e
FN<-input$n.e-input$event.e
FP<-input$event.c
TN<-input$n.c-input$event.c
data<-data.frame(TP,FN,FP,TN)
fit.roc<-reitsma(data)
summary(fit.roc)    # extract combined Sen and Spe from Here


tissue<-data[input$sampletype=="tissue",]
serum<-data[input$sampletype=="plasma",]
fit.tissue<-reitsma(tissue)
fit.serum<-reitsma(serum)
summary(fit.tissue)
summary(fit.serum)


pdf("Figure7.SROC2.pdf")
plot(fit.roc, xlim = c(0,0.7), ylim = c(0,0.75),main = "Diagnostic SROC(bivariate model) for CHD13 in NSCLC",cex.main=1,cex.lab=1,cex.axis=1,lwd=4,lty=1)
lines(sroc(fit.roc), lty = 2,lwd=4)
ROCellipse(fit.roc, lty = 2, pch = 2, d = TRUE)
points(fpr(data), sens(data), cex = 1,pch=2,lwd=4)
legend("bottomright", c("Observed data in NSCLC","SROC for NSCLC", "95%CI for NSCLC"), pch = c(2,NA,NA), lty = c(NA,2,1),cex=1.1,lwd=4)
dev.off()

# stepwise meta regression
subgroup<-c("sampletype","age","stage1","stage2","genderatio","methods","aim","Multipletarget","controldesign","Ad2Sc","primer")
collist<-na.omit(match(subgroup,colnames(input)))
metaregrlt<-list()
rlt1<-matrix(NA,length(collist),8)
for (i in 1:length(collist)){
  if(any(is.na(input[,collist[i]]))==T){
    inputomit<-input[-which(is.na(input[,collist[i]])),]
  }else{
    inputomit=input
  }
  metaomit<-metabin(event.e, n.e, event.c, n.c,data=inputomit,sm="OR", method="MH",method.tau="DL") #complab=c("NSCLC","NORMAL"),
  metaregrlt<-metareg(metaomit,~inputomit[,collist[i]]) 
  beta<-round(metaregrlt$b[2],2)
  ci.lb<-round(metaregrlt$ci.lb[2],2)
  ci.ub<-round(metaregrlt$ci.ub[2],2)
  pvalue<-(metaregrlt$pval[2])
  tau2<-round(metaregrlt$tau2,2)
  QE<-round(metaregrlt$QE,2)
  QEp<-(metaregrlt$QEp)
  R2<-(metaregrlt$R2)
  combination<-paste(colnames(input)[collist[i]],beta,"(",ci.lb,",",ci.ub,")",pvalue, sep=" ")
  print(paste(i,combination,sep=" "))
  rlt1[i,1]<-beta
  rlt1[i,2]<-ci.lb
  rlt1[i,3]<-ci.ub
  rlt1[i,4]<-pvalue  
  rlt1[i,5]<-tau2  
  rlt1[i,6]<-QE   
  rlt1[i,7]<-QEp   
  rlt1[i,8]<-R2     
}
rownames(rlt1)<-colnames(input)[collist]
colnames(rlt1)<-c("Coefficient","ci.lb","ci.ub","pvalue","tau2","QE","QE.P-value","R2")
rlt1
write.table(rlt1,file="Table.result.meta-regression.txt",sep="\t",row.names=T,col.names=NA)


colnames(input)
# conditional meta regression
metaregrlt<-metareg(meta1,~input$age+input$Ad2Sc+input$controldesign)
metaregrlt

# stepwise stratified analysis
subgroup<-c("sampletype","age","stage1","stage2","genderatio","methods","aim","Multipletarget","controldesign","Ad2Sc","primer")
collist<-match(subgroup,colnames(input))

meta1 <- metabin(event.e, n.e, event.c, n.c,data=input,studlab=author2,sm="OR", method="MH",complab=c("NSCLC","NORMAL"),method.tau="DL",byvar=input[,collist[1]])

colnames(input)

i=3
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
subdata
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="MH",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE)
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
filename
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()


i<-4
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
median<-median(input[,i],na.rm=T)
median
less<-which(subdata[,i]<65)
more<-which(subdata[,i]>=65)
subdata[less,i]<-"age<65"
subdata[more,i]<-"age>=65"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE,label.e="NSCLC",label.c="Normal")
meta2
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()

i<-5
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
median<-median(input[,i],na.rm=T)
median
less<-which(subdata[,i]<=55)
more<-which(subdata[,i]>55)
subdata[less,i]<-"stage1<49.45%"
subdata[more,i]<-"stage1>=49.45%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE,label.e="NSCLC",label.c="Normal")
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()

i<-6
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
median<-median(input[,i],na.rm=T)
median
less<-which(subdata[,i]<=median)
more<-which(subdata[,i]>median)
subdata[less,i]<-"stage <75.33%"
subdata[more,i]<-"stage(I+II)>=75.33%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()


i<-7
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
median<-median(input[,i],na.rm=T)
median
less<-which(subdata[,i]<=median)
more<-which(subdata[,i]>median)
subdata[less,i]<-"Male vs Female <69.1%"
subdata[more,i]<-"Male vs Female >=69.1%"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()


for (i in 12:15){
  valid<-! is.na(input[,i])
  subdata<-input[! is.na(input[,i]), ]
  meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author,sm="OR", method="MH",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE)
  filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
  pdf(filename)
  forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
  dev.off()
}
getwd()


i<-15
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
subdata[,i]<-as.array(as.character(subdata[,i]))
subdata[,i]
less<-which(subdata[,i]=="hom")
more<-which(subdata[,i]=="heter")
subdata[less,i]<-"Autogenous"
subdata[array(more),i]<-"Heterogeneous"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE,label.e="NSCLC",label.c="Normal")
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()



i<-17
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
median<-median(input[,i],na.rm=T)
median
less<-which(subdata[,i]<=1.9)
more<-which(subdata[,i]>1.9)
subdata[less,i]<-"Ad2Sc < 2"
subdata[more,i]<-"Ad2Sc >= 2"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE,label.e="NSCLC",label.c="Normal")
summary(meta2)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()

forest.meta(meta1)
pdf("metainf.eps")
forest(metainf(meta1, pooled="random"), comb.random=TRUE,fond=6)
dev.off()


i<-18
valid<-! is.na(input[,i])
subdata<-input[! is.na(input[,i]), ]
median<-median(input[,i],na.rm=T)
median
less<-which(subdata[,i] =="set1")
more<-which(subdata[,i] =="set2")
subdata[array(less),i]<-"set1"
subdata[array(more),i]<-"set2"
meta2 <- metabin(event.e, n.e, event.c, n.c,data=subdata,studlab=author3,sm="OR", method="I",method.tau="DL",byvar=subdata[,i],print.byvar=FALSE,label.e="NSCLC",label.c="Normal")
summary(meta2)
filename<-paste(colnames(subdata)[i],"_stratify.pdf",sep="")
pdf(filename)
forest.meta(meta2,fontsize=10,leftcols="studlab",rightcols=c("effect","ci"),print.tau2=FALSE,squaresize=1.0,plotwidth=grid::unit(4, "cm"))
dev.off()



