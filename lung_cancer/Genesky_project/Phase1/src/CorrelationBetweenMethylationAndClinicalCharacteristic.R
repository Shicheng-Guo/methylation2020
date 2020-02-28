setwd("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/result/methyClassification")

data<-read.table("Bigdata.txt",head=T,sep="\t",as.is=F,row.names=1)
head(data)
data<-RawNARemove(data)

dat<-data[seq(1,dim(data)[1],by=2),]
head(dat)
# stratified analysis

rlt<-matrixanova(dat[3:dim(dat)[2]])
write.table(rlt$pvalue,file="matrix.anova.txt",sep="\t")

x<-dat[3:dim(dat)[2]]

matrixanova<-function(x){
  # x is var matrix, 0 or 1, + or - 
  # remmber not to include sample id
  # this code can only to analysis (catogrical variable) and (catogrical or continous variable)
  
  rlt<-list()
  data<-data.frame(x)
  var<-unlist((lapply(data,function(x) length(names(table(x))))))
  catvar<-as.numeric(unlist(which(lapply(data,function(x) length(names(table(x))))<0.15*dim(data)[1])))
  nullvar1<-unlist(lapply(data[,catvar],function(x) sum(names(table(x))=="")))
  nullvar2<-var[catvar]
  nullvar2-nullvar1
  contvar<-as.numeric(unlist(which(lapply(data,function(x) length(names(table(x))))>=10)))
  
  colnames(x)[catvar]
  colnames(x)[contvar]
  x$type
  
  pvalue<-matrix(NA,nrow=dim(data)[2],ncol=dim(data)[2])
  table<-matrix(NA,nrow=dim(data)[2],ncol=dim(data)[2])
  OR<-matrix(NA,nrow=dim(data)[2],ncol=dim(data)[2])
  for(i in catvar){
    for (j in c(catvar,contvar)){
      if(j !=i ){
        tmp<-which(apply(cbind(data[,i],data[,j]),1,function(x) all(!is.na(x))))
        if(j %in% catvar & any(table(data[tmp,1],data[tmp,j])>5)){
          pvalue[i,j]<-chisq.test(data[tmp,i],data[tmp,j],simulate.p.value = TRUE)$p.value
          tmp<-table(data[tmp,i],data[tmp,j])
          table[i,j]<-paste(tmp[1,1],"/",tmp[1,2],":",tmp[2,1],"/",tmp[2,2],sep="")
          OR[i,j]<-calcOddsRatio(tmp,alpha=0.05,referencerow=1)
        }else if(j %in% contvar){
          pvalue[i,j]<-try(oneway.test(as.numeric(data[tmp,j]) ~ data[tmp,i],na.action="na.omit")$p.value)
          print (c(i,j))
        }
      }
    }
  }

  colnames(pvalue)=colnames(data)
  rownames(pvalue)=colnames(data)
  colnames(OR)=colnames(data)
  rownames(OR)=colnames(data)
  colnames(table)=colnames(data)
  rownames(table)=colnames(data)
  rlt$pvalue=pvalue
  rlt$table=table
  rlt$OR=OR
  return(rlt)
}



dat1<-data[which(data$Smoking>0),c(3,4,5,10,15,16,17,18,23,9)]
dat2<-data[which(data$Smoking==0),c(3,4,5,10,15,16,17,18,23,9)]

dat3<-data[which(data$type=="ad"),c(3,4,5,10,15,16,17,18,23,9)]
dat4<-data[which(data$type=="sq"),c(3,4,5,10,15,16,17,18,23,9)]

dat5<-data[which(data$differention=="low"),c(3,4,5,10,15,16,17,18,23,9)]
dat6<-data[c(which(data$differention=="medium"),which(data$differention=="high")),c(3,4,5,10,15,16,17,18,23,9)]

early<-c(which(data$TNM=="IA"),which(data$TNM=="IB"),which(data$TNM=="IIA"),which(data$TNM=="IIB"))
late<-c(which(data$TNM=="IIIA"),which(data$TNM=="IIIB"),which(data$TNM=="IV"))
dat7<-data[early,c(3,4,5,10,15,16,17,18,23,9)]
dat8<-data[late,c(3,4,5,10,15,16,17,18,23,9)]

head(data)

