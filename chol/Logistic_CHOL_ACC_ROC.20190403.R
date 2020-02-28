install.packages("pROC")
install.packages("ggplot2")

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

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}

options(digits = 2)
index2type<-function(index){
  sampletype=ifelse(as.numeric(index) %% 2,"Case","Normal") # for chol project, odds is case while even is control
}
data2summary <- function(data, varname, groupnames){
  # require(plyr)
  # c(mean(x)-2*sem,mean(x)+2*sem)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem=sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])),
      iqr=as.numeric(quantile(x[[col]],na.rm=T)[4]-quantile(x[[col]],na.rm=T)[2]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Table2Generator = function(methdata){
  seq.case = which(methdata[,1] ==1)
  seq.control = which(methdata[,1] == 0)
  #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methdata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case],na.rm=T)$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr")
  #Logistic regression analysis
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:(ncol(methdata)-1 )){
    temp = methdata[,c(1,i+1 )]
    temp[,1] = ifelse(temp[,1] ==1,1,0)
    temp[,1] = as.factor(temp[,1])
    glm.fit  = glm(temp[,1] ~ temp[,2], data = temp, family = "binomial")
    OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
    Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
    CI.upper[i]=log(exp(confint(glm.fit)[2,2]),base = 10)
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
    #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit)
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]]
    #Find the best Sens and Spec
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = logistic.rocdata[seq.max,1]
    Spec[i] = logistic.rocdata[seq.max,2]
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}

rlt<-combineAUC(methdata,recombination=c(i))
recombination=i

combineAUC<-function(methdata,recombination="."){
  Table<-list()
  temp <- methdata[,grepl(paste(recombination, collapse="|"), colnames(methdata))]
  genesymbol= unlist(lapply(colnames(temp), function(x) strsplit(as.character(x),"_")[[1]][1]))
  temp<-t(apply(temp,1,function(x) tapply(x, genesymbol,function(x) mean(x,na.rm=T))))
  head(temp)
  if(nrow(temp)==1){
    temp<-t(temp)
    colnames(temp)<-unique(genesymbol)
  }
  head(temp)
  phen=as.numeric(rownames(methdata)) %%2   # 0 control, 1 case
  newinput=data.frame(phen,temp)
  head(newinput)
  glm.fit  = glm(phen~ ., data = newinput, family = "binomial")
  summary(glm.fit)
  logOR = log(exp(summary(glm.fit)$coefficients[,1]),base = 10)
  Logistic.P = summary(glm.fit)$coefficients[,4]
  CI.upper=log(exp(confint(glm.fit)[,2]),base = 10)
  CI.lower = log(exp(confint(glm.fit)[,1]),base = 10)
  Mean<-tapply(newinput[,2],newinput[,1],function(x) mean(x,na.rm=T))
  SD<-tapply(newinput[,2],newinput[,1],function(x) sd(x,na.rm=T))
  #Do the analysis of the sens, spec, and AUC
  predicted.value = predict(glm.fit)
  
  pred <- predict(glm.fit,newinput,type="response")
  real <- newinput$phen
  plot.roc(real,pred, col = 3, main="ROC Validation set",percent = TRUE, print.auc = TRUE)
  
  predicted.data  = data.frame(Type=na.omit(newinput)[,1], predicted.value)
  logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
  logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
  AUC = logistic.rocobj$auc[[1]]
  #Find the best Sens and Spec
  logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
  seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
  Sens = logistic.rocdata[seq.max,1]
  Spec = logistic.rocdata[seq.max,2]
  Table$matrix = data.frame(logOR, CI.upper, CI.lower, Logistic.P)
  Table$model=c(MFO=Mean[1], MFC=Mean[2],SD0=SD[1],SD1=SD[2],logOR=logOR[2],Pval=Logistic.P[2],CI_upper=CI.upper[2],CI_lower=CI.lower[2],Sen=Sens,Spec=Spec,AUC=AUC)
  Table$roc=logistic.rocobj
  return(Table)
}


bestcombineAUC<-function(methdata,recombination="."){
  Table<-list()
  temp <- methdata[,grepl(paste(recombination, collapse="|"), colnames(methdata))]
  genesymbol= unlist(lapply(colnames(temp), function(x) strsplit(as.character(x),"_")[[1]][1]))
  temp<-t(apply(temp,1,function(x) tapply(x, genesymbol,function(x) mean(x,na.rm=T))))
  head(temp)
  if(nrow(temp)==1){
    temp<-t(temp)
    colnames(temp)<-unique(genesymbol)
  }
  head(temp)
  phen=rep(0,nrow(temp))  
  phen[grep("T",rownames(temp))]<-1
  
  temp=na.omit(data.frame(phen,temp))
  
  glm.null <- glm(phen ~ 1, data = temp,family = "binomial")
  glm.fit  = glm(phen~ ., data = temp, family = "binomial")
  step_model <- step(glm.null, scope = list(lower = glm.null, upper = glm.fit), direction = "forward")
  
  summary(step_model)
  pred <- predict(step_model,temp[,2:ncol(temp)],type="response")
  real <- temp$phen
  plot.roc(real,pred, col = 3, main="ROC Validation set",percent = TRUE, print.auc = TRUE)
  
  summary(step_model)
  logOR = log(exp(summary(step_model)$coefficients[,1]),base = 10)
  Logistic.P = summary(step_model)$coefficients[,4]
  CI.upper=log(exp(confint(step_model)[,2]),base = 10)
  CI.lower = log(exp(confint(step_model)[,1]),base = 10)
  #Do the analysis of the sens, spec, and AUC
  predicted.value = predict(step_model)
  predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
  logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
  logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
  AUC = logistic.rocobj$auc[[1]]
  #Find the best Sens and Spec
  logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
  seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
  Sens = logistic.rocdata[seq.max,1]
  Spec = logistic.rocdata[seq.max,2]
  Table$matrix = data.frame(logOR, CI.upper, CI.lower, Logistic.P)
  Table$model=c(Sen=Sens,Spec=Spec,AUC=AUC)
  Table$roc=logistic.rocobj
  return(Table)
}

options(digits = 2)
index2type<-function(index){
  sampletype=ifelse(as.numeric(index) %% 2,"Case","Normal") # for chol project, odds is case while even is control
}

data2summary <- function(data, varname, groupnames){
  # require(plyr)
  # c(mean(x)-2*sem,mean(x)+2*sem)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem=sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])),
      iqr=as.numeric(quantile(x[[col]],na.rm=T)[4]-quantile(x[[col]],na.rm=T)[2]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
#############################################################################################################
#############################################################################################################
#############################################################################################################

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-1/profile")
library(pROC)
library(readxl)
library(ggplot2)

data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))

methdata<-data.matrix(data[,c(7:213)])
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
phen=as.numeric(colnames(methdata)) %%2   # 0 control, 1 case
methdata=data.frame(phen,t(methdata))
methdata[1:5,1:5]

methdata=ColNARemove(methdata)
methdata[1:5,1:5]
genesymbol= unlist(lapply(colnames(methdata)[2:ncol(methdata)], function(x) strsplit(as.character(x),"_")[[1]][1]))
head(genesymbol)

rlt<-Table2Generator(methdata)
genesymbol= unlist(lapply(rownames(rlt), function(x) strsplit(as.character(x),"_")[[1]][1]))
rlt<-data.frame(genesymbol,rlt)

write.table(rlt,file="Phase2.CHOL.Sen.Spe.Acc.txt",col.names=NA,row.names = T,quote=F,sep="\t")
xueqing<-rlt[rlt$genesymbol %in% c("SOX11","LINE.1","SHOX2"),]
write.table(xueqing,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.txt",col.names=NA,row.names = T,quote=F,sep="\t")

rlt1<-combineAUC(methdata,c("LINE.1","SHOX2"))
rlt2<-combineAUC(methdata,c("LINE.1","SOX11"))
rlt3<-combineAUC(methdata,c("SHOX2","SOX11"))
rlt4<-combineAUC(methdata,c("LINE.1","SHOX2","SOX11"))

write.table(rlt1$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt1.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt1$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt1.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt2$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt2.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt2$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt2.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt3$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt3.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt3$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt3.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt4$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt4.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt4$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt4.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

rlt6<-combineAUC(methdata,recombination=c("LINC00466","LncRNA"))
rlt7<-combineAUC(methdata,recombination=c("LINC00466","MIR129"))
rlt8<-combineAUC(methdata,recombination=c("MIR129","LncRNA"))
rlt9<-combineAUC(methdata,recombination=c("LINC00466","MIR129","LncRNA"))

write.table(rlt6$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt6.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt6$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt6.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt7$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt7.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt7$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt7.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt8$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt8.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt9$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt8.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt9$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt9.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt9$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt9.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

plot(rlt6$roc,col=3)
lines(rlt7$roc,col=4)
lines(rlt8$roc,col=5)
lines(rlt9$roc,col=2)
legend("bottomright",cex=0.9,legend=c("LINC00466+TCONS_00021941","LINC00466+MIR129-2","MIR129-2+TCONS_00021941","LINC00466+MIR129-2+TCONS_00021941"),col=c(3,4,5,2),lty=1,bty="n")

rlt10<-combineAUC(methdata,recombination=c("LINC00466"))
rlt11<-combineAUC(methdata,recombination=c("SOX11"))
rlt12<-combineAUC(methdata,recombination=c("LINE.1"))
rlt13<-combineAUC(methdata,recombination=c("SHOX2"))

pdf("SOX11.LINC00466.pdf")
par(mfrow=c(2,2))
plot(rlt10$roc,col=3,main="LINC00466")
text(x=0.25,y=0.25,paste("AUC=",round(rlt10$model[3],3),sep=""))
plot(rlt11$roc,col=3,main="SOX11")
text(x=0.25,y=0.25,paste("AUC=",round(rlt11$model[3],3),sep=""))
plot(rlt12$roc,col=3,main="LINE-1")
text(x=0.25,y=0.25,paste("AUC=",round(rlt12$model[3],3),sep=""))
plot(rlt13$roc,col=3,main="SHOX2")
text(x=0.25,y=0.25,paste("AUC=",round(rlt13$model[3],3),sep=""))
dev.off()

AUC<-c()
GENE<-unlist(lapply(unique(genesymbol),function(x) unlist(strsplit(x,"[-]"))[1]))
head(GENE)
for(i in GENE){
  print(i)
  rlt<-combineAUC(methdata,recombination=c(i))
  AUC<-rbind(AUC,rlt$model)
}
rownames(AUC)<-unique(genesymbol)
AUC
write.table(AUC,file="Phase2.CHOL.Sen.Spe.Average.AUC.txt",col.names=NA,row.names = T,quote=F,sep="\t")
rlt14<-combineAUC(methdata,recombination=c("SHOX2","CDO1","SLC6A2","TRH","AKAP12"))
rlt14
genesymbol= unlist(lapply(colnames(methdata)[2:ncol(methdata)], function(x) strsplit(as.character(x),"_")[[1]][1]))
rlt15<-combineAUC(methdata,recombination=unique(genesymbol))
rlt15
rlt14<-combineAUC(methdata,recombination=c("CYP2A13","ZNF415","NTM","LINE.1"))
rlt14
######################################################################################################################
##################### 2019-02-22 (HCC9810 and RBE with 5-AZA treatment) ##############################################
######################################################################################################################
install.packages("ggplot2")
install.packages("colorspace")
library("colorspace")
library("plyr")
library("ggplot2")

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-1/profile")

data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
colnames(data)
methdata<-data.matrix(data[,c(7:ncol(data))])
head(methdata)
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rowname)
# even number is control while odd number is case
sampletype=as.numeric(colnames(methdata)) %%2   # 0 control, 1 case
mydata<-list()
mydata$methdata=methdata
mydata$sampletype<-sampletype
mydata$genomeposition<-data$GenomePosition
mydata$genesymbol<-genesymbol
# wide format to long format
library(tidyr)
newdata<-data.frame(GeneSymbol=genesymbol,GenomePosition=mydata$genomeposition,mydata$methdata,check.names=F)
keycol <- "index"
valuecol <- "mf"
gathercols <- colnames(newdata)[3:ncol(newdata)]
newdataLong<-gather_(newdata, keycol, valuecol, gathercols)
head(newdataLong)
newdataLong$type=as.numeric(newdataLong$index) %%2   # 0 control, 1 case
head(newdataLong)
input <- data2summary(newdataLong, varname="mf",groupnames=c("GeneSymbol","GenomePosition","type")) 
head(input)

for(i in unique(input$GeneSymbol)){
  options(scipen=10000)
  inputsubset<-subset(input,GeneSymbol==i)
  p<- ggplot(inputsubset, aes(x=GenomePosition, y=mf, group=type, color=as.factor(type))) + 
    geom_line(size=1.5) +
    geom_point()+
    geom_errorbar(aes(ymin=mf-2*sem, ymax=mf+2*sem), width=.2,size=1.5,position=position_dodge(0.05))
  p<-p+labs(title=paste("methylation profile of",i,sep=" "), x="Position (bp)", y = "Methylation level (beta)")+
    theme_classic() +
    scale_color_manual(values=c('red',"blue"))+
    theme(axis.line = element_line(color = 'black'))+
    theme(axis.text.x=element_text(size=15,vjust=0.9,face="bold",hjust=1,angle=45))+
    theme(axis.title.x=element_text(size=25,vjust=-0.1))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.y=element_text(size=22,vjust=0.3)) #draws x and y axis line
    ggsave(paste(i,".profile.pdf",sep=""))
}
# Prepare ggplot matrix for each gene/region
# Split the positive control datasets
library("beeswarm")
for(i in unique(newdataLong$GeneSymbol)){
  print(i)
  pdf(paste(i,".beeswarm.pdf",sep=""))
  options(scipen=10000)
  newdataLongsubset<-subset(newdataLong,GeneSymbol==i)
  plot.data<-data.frame(mf=newdataLongsubset$mf,type=newdataLongsubset$type)
  try(beeswarm(mf ~ type, data = plot.data, pch = 16, xlab = "", ylab = "Methylation level (beta)",labels = unique(newdataLong$type), col= 1:length(unique(newdataLong$type))))
  legend("topright", legend = unique(newdataLong$type),title = "Sample", pch = 16, col = 1:length(unique(newdataLong$type)))
  dev.off()
}
###############################################################################################################################
###########################                                                        ############################################
###############################################################################################################################
