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
library("tidyr")

options(digits = 2)
index2type<-function(index){
  sampletype=ifelse(substr(index,1,1)=="T","Case","Normal") # for chol project, odds is case while even is control
  sampletype
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

install.packages("xlsx")
library("xlsx")

setwd("C:\\Users\\shg047\\Dropbox\\Realworld")
list.files()
for(i in 13:20){
  data= read_excel("RA_Realworld_EDC_20180306.xlsx",sheet =i)
  write.xlsx(data,file=paste("RA_Realworld_EDC_20180306_Sheet",i,".xlsx",sep=""), sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)
  print(i)
}
setwd("/home/shg047/Dropbox/Realworld/")
list.files()
data= read_excel("RA_Realworld_EDC_20180306_Sheet15.xlsx",sheet =1,skip=1)
head(data)
data<-as.data.frame(data)
head(data)

data= read_excel("RA_Realworld_EDC_20180306_Sheet9.xlsx",sheet =1,skip=1)
head(data)
data<-as.data.frame(data)
head(data)

