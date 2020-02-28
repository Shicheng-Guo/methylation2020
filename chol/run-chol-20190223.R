library("ggplot2")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/result")

data1= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-1/methylation.xlsx",sheet = 2)
data2= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-2/methylation.xlsx",sheet = 2)
data1= as.data.frame(data1)
data2= as.data.frame(data2)
data1<-data1[,1:178]
data2<-data2[,1:48]
rowname1<-apply(data.frame(data1$Target,as.character(data1$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
rowname2<-apply(data.frame(data2$Target,as.character(data2$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
rowname1==rowname2
rowname<-rowname1


methdata1<-data.matrix(data1[,c(7:176)])
methdata2<-data.matrix(data2[,c(7:48)])
colnames(methdata1)
colnames(methdata2)<-as.numeric(colnames(methdata2))+200
rownames(methdata1)<-rowname
rownames(methdata2)<-rowname
GenomePosition<-c(data1$GenomePosition)

methdata<-data.frame(methdata1,methdata2,check.names=F)
dim(methdata)
colnames(methdata)
methdata[1:5,1:5]

genesymbol= unlist(lapply(data1$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rowname)
phen=as.numeric(gsub("X","",colnames(methdata))) %%2   # 0 control, 1 case
head(phen)
head(colnames(methdata))

methdata=data.frame(phen,t(methdata))
methdata[1:5,1:5]

Table2Generator = function(methydata){
  seq.case = which(methydata[,1] ==1)
  seq.control = which(methydata[,1] == 0)
  # Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methydata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methydata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methydata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case])$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr")
  
  #Logistic regression analysis
  library(pROC)
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:(dim(methydata)[2] -1 )){
    temp = methydata[,c(1,i+1 )]
    glm.fit  = glm(temp[,1] ~ temp[,2], data = temp, family = "binomial")
    OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
    Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
    CI.upper[i] = log(exp(confint(glm.fit)[2,2]),base = 10)
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
    #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit,type=c("response"))
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]]
    #Find the best Sens and Spec
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = logistic.rocdata[seq.max,1]
    Spec[i] = logistic.rocdata[seq.max,2]
    print(paste(i,colnames(methydata)[i]))
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}

Table<-Table2Generator(methdata)
write.table(Table,file="chol.rlt.txt",sep="\t",quote=F,col.names = NA,row.names =T)