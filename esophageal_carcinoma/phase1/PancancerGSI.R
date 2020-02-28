########################################################################################
###   Title: Tissue specific hyper-methylation biomarkers
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 11/9/2015
###   Strategy: 1) Plasma based Cancer biomaker require tissue specific hyper-methylation. 
###   Strategy: 2) Methylation Biomarker should be low and non-methylated in PBMC samples
########################################################################################

###########################################################################################
########################### First Step: Load R Function ###################################
###########################################################################################

gsi<-function(data){
  gsi<-c()
  gmaxgroup<-c()
  cpg<-c()
  for(i in 1:nrow(data)){
    newdata<-t(na.omit(t(data[i,])))
    group=names(table(colnames(newdata)))
    index=colnames(newdata)
    newdata<-as.numeric(newdata)
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(newdata),index,mean)))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(newdata[which(index==group[j])]))/10^(mean(newdata[which(index==gmax)])))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    gsi<-c(gsi,gsit)
    # print(c(i,rownames(data)[i],gmax,gsit))
  }
  rlt=data.frame(CpG=rownames(data),group=gmaxgroup,GSI=gsi)
  return(rlt)
}

RawNARemove<-function(data,missratio=1){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
}   

# ========================================================================================
args <- commandArgs(trailingOnly = TRUE)
print(args)
data<-read.table(args[1],sep="\t",colClasses=c("character",rep("numeric",6440)),nrow=2500,row.names=1)
colname<-read.table(args[2],sep="\t",row.names=1,colClasses=c("character",rep("character",6440)),nrow=1)
saminfo<-read.table(args[3],sep="\t")
group<-paste(saminfo[match(colname,saminfo[,1]),3],"_",saminfo[match(colname,saminfo[,1]),4],sep="")
colnames(data)=group
data=RawNARemove(data)
rlt<-gsi(data)
write.table(rlt,file=args[4],sep="\t",col.names=F,row.names=F,quote=F)

