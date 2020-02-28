source("https://bioconductor.org/biocLite.R")
if (!require("ggsci")) biocLite("ggsci")
if (!require("readxl")) biocLite("readxl")
if (!require("readr")) biocLite("readr")
if (!require("org.Hs.eg.db")) biocLite("org.Hs.eg.db")
if (!require("GO.db")) biocLite("GO.db")
if (!require("ggbeeswarm")) biocLite("ggbeeswarm")
if (!require("gridExtra")) biocLite("gridExtra")
if (!require("biomaRt")) biocLite("biomaRt")
if (!require("knitr")) biocLite("knitr")
if (!require("ggplot2")) biocLite("ggplot2")
if (!require("ggthemes")) biocLite("ggthemes")
if (!require("reshape2")) biocLite("reshape2")
if (!require("GOstats")) biocLite("GOstats")
if (!require("pROC")) install.packages("pROC_1.13.0.tar.gz")
install.packages("rmarkdown")


library("knitr")
library("ggsci")
library("ggplot2")
library("ggthemes")
library("readxl")
library("reshape2")
library("readr")
library("dplyr")
library("org.Hs.eg.db")
library("GO.db")
library("GOstats")
library("ggbeeswarm")
library("gridExtra")
library("biomaRt")
library("reshape2")
library("plyr")
library("pROC")

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

#####################################################################################################
############### 2019-02-22 (HCC9810 and RBE with 5-AZA treatment) ##################################
#####################################################################################################
library("ggplot2")
# setwd("/media/Home_Raid1/shg047/NAS1/Minghua2017/chol")
setwd("/home/guosa/hpc/methylation/chol/final")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))

methdata<-data.matrix(data[,c(181:184)])
head(methdata)
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rowname)
data[1:5,1:10]
# even number is control while odd number is case
sampletype=colnames(data)[181:184]
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
gathercols <- colnames(methdata)
newdataLong<-gather_(newdata, keycol, valuecol, gathercols)
newdataLong$type=colnames(data)[181:184]
input <- data2summary(newdataLong, varname="mf",groupnames=c("GeneSymbol","GenomePosition","type")) 
head(input)

for(i in unique(input$GeneSymbol)){
  pdf(paste(i,".profile.aza.pdf",sep=""))
  options(scipen=10000)
  inputsubset<-subset(input,GeneSymbol==i)
  p<- ggplot(inputsubset, aes(x=GenomePosition, y=mf, group=type, color=type)) + 
    geom_line(size=1.5) +
    geom_point()+
    geom_errorbar(aes(ymin=mf-2*sem, ymax=mf+2*sem), width=.2,size=1.5,position=position_dodge(0.05))
  p<-p+labs(title=paste("methylation profile of",i,sep=" "), x="Position (bp)", y = "Methylation level (beta)")+
    theme_classic() +
    scale_color_manual(values=c('red','blue',"green","yellow"))+
    theme(axis.line = element_line(color = 'black'))+
    theme(axis.text.x=element_text(size=15,vjust=0.9,face="bold",hjust=1,angle=45))+
    theme(axis.title.x=element_text(size=25,vjust=-0.1))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.y=element_text(size=22,vjust=0.3)) #draws x and y axis line
  print(p)
  dev.off()
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

#####################################################################################################
##########################                                   ########################################
#####################################################################################################
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/phase2")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
head(rowname)
data[1:10,1:10]

methdata<-data.matrix(data[,-c(1:6)])
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rowname)
data[1:10,1:10]
# even number is control while odd number is case
sampletype=ifelse(as.numeric(colnames(methdata)) %% 2,"Case","Normal")
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
gathercols <- colnames(methdata)
newdataLong<-gather_(newdata, keycol, valuecol, gathercols)
newdataLong$type=index2type(newdataLong$index)
input <- data2summary(newdataLong, varname="mf",groupnames=c("GeneSymbol","GenomePosition","type")) 

for(i in unique(input$GeneSymbol)){
  pdf(paste(i,".profile.pdf",sep=""))
  options(scipen=10000)
  inputsubset<-subset(input,GeneSymbol==i)
  p<- ggplot(inputsubset, aes(x=GenomePosition, y=mf, group=type, color=type)) + 
    geom_line(size=1.5) +
    geom_point()+
    geom_errorbar(aes(ymin=mf-2*sem, ymax=mf+2*sem), width=.2,size=1.5,position=position_dodge(0.05))
  p<-p+labs(title=paste("methylation profile of",i,sep=" "), x="Position (bp)", y = "Methylation level (beta)")+
    theme_classic() +
    scale_color_manual(values=c('red','blue'))+
    theme(axis.line = element_line(color = 'black'))+
    theme(axis.text.x=element_text(size=15,vjust=0.9,face="bold",hjust=1,angle=45))+
    theme(axis.title.x=element_text(size=25,vjust=-0.1))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.y=element_text(size=22,vjust=0.3)) #draws x and y axis line
  print(p)
  dev.off()
}
# Prepare ggplot matrix for each gene/region
# Split the positive control datasets
library("beeswarm")
for(i in unique(newdataLong$GeneSymbol)){
  pdf(paste(i,".beeswarm.pdf",sep=""))
  options(scipen=10000)
  newdataLongsubset<-subset(newdataLong,GeneSymbol==i)
  plot.data<-data.frame(mf=newdataLongsubset$mf,type=newdataLongsubset$type)
  beeswarm(mf ~ type, data = plot.data, pch = 16, xlab = "", ylab = "Methylation level (beta)",labels = c("Cancer", "Normal"), col=2:3)
  legend("topright", legend = c("Cancer", "Normal"),title = "Sample", pch = 16, col = 2:3)
  dev.off()
  print(i)
}
###########################################################################################################
########################                                                    ###############################
###########################################################################################################
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/phase2")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
methdata<-data.matrix(data[,c(7:176)])
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rowname)
phen=as.numeric(colnames(methdata)) %%2   # 0 control, 1 case
methdata=data.frame(phen,t(methdata))
methdata[1:5,1:5]

Table2Generator = function(methydata){
  seq.case = which(methydata[,1] ==0)
  seq.control = which(methydata[,1] == 1)
  #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methydata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methydata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methydata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case])$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr")
  #Logistic regression analysis
  library("pROC")
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

methdata[1:5,1:5]
GeneSymbol<-unique(sapply(colnames(methdata)[2:ncol(methdata)],function(x) unlist(strsplit(x,split="[_]"))[[1]]))
GeneSymbol<-GeneSymbol[-grep("cg",GeneSymbol)]
OR<-c()
Logistic.P<-c()
CI.upper<-c()
CI.lower<-c()
Group<-c()
for(i in 1:length(GeneSymbol)){
  data<-methdata[,grep(GeneSymbol[i],colnames(methdata))]
  newdata<-apply(data,1,function(x) mean(x,na.rm=T))
  phen=as.numeric(rownames(data)) %%2   # 0 control, 1 case
  glm.fit  = glm(phen ~ newdata, family = "binomial")
  OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
  Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
  CI.upper[i] = log(exp(confint(glm.fit)[2,2]),base = 10)
  CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
  Group<-c(Group,GeneSymbol)
}
Table = data.frame(OR,LowerLimit=CI.lower,UpperLimit=CI.upper,Group)

p = ggplot(data=Table,
           aes(x = Group,y = OR, ymin = LowerLimit, ymax = UpperLimit))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab('Group')+ ylab("Odds Ratio (95% Confidence Interval)")+
  coord_flip() + geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.5,cex=1)+ 
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(size=10,face="bold"),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  theme(legend.position="none")
p
