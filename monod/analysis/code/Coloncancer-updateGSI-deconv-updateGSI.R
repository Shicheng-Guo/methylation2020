Coloncancerdatapre<-function(){
  rlt<-list()
  setwd("/home/shg047/oasis/monod/hapinfo/")
  saminfo<-read.table("/home/shg047/oasis/monod/saminfo.txt",sep="\t")
  data<-read.table("/home/shg047/oasis/monod/hapinfo/monod.mhl.may5.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
  data<-data[,-grep("CTT-|PC-P|PC-T",colnames(data))]
  colnames(data)[grep("STL",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("STL",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("WB",colnames(data))]<-"WBC"
  colnames(data)[grep("N37",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("N37",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("methylC",colnames(data))]<-"H1"
  colnames(data)[grep("6-T|Colon_Tumor_Primary",colnames(data))]<-"CCT"
  colnames(data)[grep("7-T-",colnames(data))]<-"LCT"
  colnames(data)[grep("6-P-",colnames(data))]<-"CCP"
  colnames(data)[grep("7-P-",colnames(data))]<-"LCP"
  colnames(data)[grep("NC-P-",colnames(data))]<-"NCP"
  data1<-data[,grep("Brain|Stomach|Lung|Heart|CCT|Colon|WBC|Liver|Esophagus|Kidney|CCP",colnames(data))]
  data1_ref<-data1[,-(grep("CCP",colnames(data1)))]
  base=unlist(lapply(strsplit(colnames(data1_ref),"[.]"),function(x) x[[1]]))
  colnames(data1_ref)<-base
  gsi<-data.frame(GSIMethDecon(data1_ref))
  group<-as.character(unique(gsi$group))
  signatures<-AverageWithGroup(data1_ref)  # QR
  rlt$GSI<-gsi
  rlt$data<-data1
  rlt$ref<-data1_ref
  rlt$signature<-signatures
  return(rlt)
}

GSIMethDecon<-function(data){
  data<-data.matrix(data)
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-gmaxgroup<-avebase<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(data[i,][which(index==group[j])],na.rm=T))/10^(mean(data[i,][which(index==gmax)],,na.rm=T)))/(length(group)-1)
      gsit<-gsit+tmp
    }
    ave<-tapply(data[i,], index, function(x) mean(x,na.rm=T))
    gmaxgroup<-c(gmaxgroup,gmax)
    gsi<-c(gsi,gsit)
    avebase<-rbind(avebase,ave)
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi,AVE=avebase)
  rlt<-rlt[-which(rlt[,2]!="WBC" & rlt[,grep("WBC",colnames(rlt))]>0.01),]
  return(rlt)
}


AverageWithGroup<-function(data){
  base=unlist(lapply(strsplit(colnames(data),"[.]"),function(x) x[[1]]))
  matrix=apply(data,1,function(x) tapply(x,base,function(x) mean(x,na.rm=T)))
  matrix<-t(matrix)
  rownames(matrix)=rownames(data)
  matrix<-matrix[!rowSums(!is.finite(matrix)),]
  return(matrix)
}

RandomSamplingMean<-function(data,number=round(ncol(data)/2)){
  rlt<-c()
  for(i in 1:ncol(data)){
    VirtualSample<-rowMeans(data[,sample(1:ncol(data),number)],na.rm=T)  
    rlt<-cbind(rlt,VirtualSample)
  }
  rlt
}

TopGSIByCategory<-function(gsi,top=20){
  GSIRlt<-c()
  rank<-c(rep(top,length(group)))
  for (i in 1:length(group)){
    subset=gsi[which(gsi$group==group[i]),]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
  }
  return(GSIRlt)
}


DataExtractionTrim<-function(data){
  data<-data[,-grep("CTT-|PC-P|PC-T",colnames(data))]
  colnames(data)[grep("STL",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("STL",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("WB",colnames(data))]<-"WBC"
  colnames(data)[grep("N37",colnames(data))]<-as.character(saminfo[match(colnames(data)[grep("N37",colnames(data))],saminfo[,1]),2])
  colnames(data)[grep("methylC",colnames(data))]<-"H1"
  colnames(data)[grep("6-T|Colon_Tumor_Primary",colnames(data))]<-"CCT"
  colnames(data)[grep("7-T-",colnames(data))]<-"LCT"
  colnames(data)[grep("6-P-",colnames(data))]<-"CCP"
  colnames(data)[grep("7-P-",colnames(data))]<-"LCP"
  colnames(data)[grep("NC-P-",colnames(data))]<-"NCP"
  data1<-data[,grep("Brain|Stomach|Lung|Heart|Colon|CCT|WBC|Liver|Esophagus|Kidney|Intestine|CCP",colnames(data))]
  return(data1)
}

FigurePrepareSimulation<-function(){
  library("ggplot2")
  Fig<-data.frame(seq(0,1,by=0.05),Rlt$out.all[,grep("CCT",colnames(Rlt$out.all))],Rlt$out.all[,grep("WB",colnames(Rlt$out.all))])
  colnames(Fig)<-c("Simulated","CCT","WB")
  Fig<-100*Fig
  c <- ggplot(Fig, aes(Simulated,CCT))
  c <- c + xlab("Simulated contribution (%)")+ylab("Predicted contribution (%)")
  c <- c + stat_smooth(se = TRUE,n=10,size=1.5) + geom_point(size=3)
  c <- c + stat_smooth(se = TRUE,n=10,size=1.5) + geom_point(size=3)
  c <- c + theme_bw()
  c <- c + theme(axis.text=element_text(size=10), axis.title=element_text(size=14,face="bold"))
  ggsave(c,file="coloncancer-deconvolution.simultaion-2.pdf")
  save.image(file="Colon-Deconvolution.RData")
}


ci95<-function(x){
  x<-data.frame(x)
  out<-c()
  for(i in 1:ncol(x)){
    error <- qt(0.975,df=length(x[,i])-1)*sd(x[,i])/sqrt(length(x[,i]))
    m<-round(mean(x[,i]),4)
    d<-round(mean(x[,i])-error,4)
    u<-round(mean(x[,i])+error,4)
    rlt<-paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
    out<-cbind(out,rlt)
  }
  return(out)
}


library("DeconRNASeq")
library("ggplot2")

Colon<-Coloncancerdatapre()
data1=Colon$data
data1_ref=Colon$ref
signatures<-AverageWithGroup(data1_ref)  
topgsi<-TopGSIByCategory(Colon$GSI,top=50)
VirtualMatrix<-data.frame(RandomSamplingMean(data1[,grep("CCP",colnames(data1))]))
VirtualMatrix<-VirtualMatrix[-which(! is.finite(rowSums(VirtualMatrix))),]
head(VirtualMatrix)
# Deconvolution-QR
DeconData<-data.frame(VirtualMatrix,signatures[match(rownames(VirtualMatrix),rownames(signatures)),])
DeconData.tmp<-DeconData[na.omit(match(topgsi[,1],rownames(DeconData))),]
VirtualMatrix=data.frame(DeconData.tmp[,grep("VirtualSample",colnames(DeconData))])
Signatures=data.frame(DeconData.tmp[,-grep("VirtualSample",colnames(DeconData))])
library(car)
NewVirtualMatrix=logit(VirtualMatrix)
NewSignatures=logit(Signatures)
head(NewVirtualMatrix)
head(NewSignatures)
Rlt<-DeconRNASeq(NewVirtualMatrix,NewSignatures,checksig=FALSE,known.prop = F, use.scale = TRUE, fig = TRUE)
write.table(Rlt$out.all,file="Realdata-Colon.deconvolution.txt",sep="\t",quote=F,col.names=NA,row.names=T)
save.image(file="colon.cancer.updateGSI0.01.RData")


