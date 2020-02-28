library("pROC")
library("ggplot2")

Table2Generator = function(methdata){
  seq.case = which(methdata[,1] ==1)
  seq.control = which(methdata[,1] == 11)
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
    print(i)
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}

id2phen4<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
}

id2phen3<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*"))
}

id2bin<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
  as.numeric(lapply(strsplit(filename,"-"),function(x) x[4]))
}

id2pid<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"edu_...."))
  unlist(lapply(filename,function(x) unlist(strsplit(x,"[_]"))[2]))
}

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

Table2Generator = function(methdata){
  seq.case = which(methdata[,1] ==1)
  seq.control = which(methdata[,1] == 11) #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methdata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case],na.rm=T)$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr") #Logistic regression analysis
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
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10) #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit)
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]] #Find the best Sens and Spec
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
  # require(plyr);  # c(mean(x)-2*sem,mean(x)+2*sem)
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
  phen=rep(0,nrow(temp))  
  phen[grep("T",rownames(temp))]<-1
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




setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
methdata[1:5,1:5]

methdata<-methdata[,grep("LGG",colnames(methdata))]
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

GBM<-grep("GBM",colnames(input))
newinput<-input[,GBM]
newphen<-phen[GBM,]
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

HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  hr<-summary(coxph(Surv(month,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
rownames(HR)<-rownames(newdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
rlt<-data.frame(HR,map[match(rownames(HR),map[,4]),])
write.table(rlt,file="~/hpc/methylation/TCGA_HM450_GBM_Suvival_HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)




setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
methdata[1:5,1:5]
methdata<-methdata[,grep("LGG",colnames(methdata))]
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
HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  hr<-summary(coxph(Surv(month,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
rownames(HR)<-rownames(newdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
rlt<-data.frame(HR,map[match(rownames(HR),map[,4]),])
write.table(rlt,file="~/hpc/methylation/TCGA_HM450_LGG_Suvival_HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)





