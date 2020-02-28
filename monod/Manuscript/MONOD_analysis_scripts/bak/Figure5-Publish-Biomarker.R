  gsi<-function(data){
    group=names(table(colnames(data)))
    index=colnames(data)
    GSI<-c()
    gmaxgroup<-c()
    for(i in 1:nrow(data)){
      gsit<-0
      gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
      if(length(gmax)<1){print(data[i,])}
      for(j in 1:length(group)){
        tmp<-(1-10^(mean(na.omit(as.numeric(data[i,which(index==group[j])])),na.rm=T))/10^(mean(na.omit(as.numeric(data[i,which(index==gmax)])))))/(length(group)-1)
        gsit<-c(gsit,tmp)
      }
      print(i)
      gmaxgroup<-c(gmaxgroup,gmax)
      GSI<-c(GSI,sum(gsit,na.rm=T))
    }
    rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=GSI)
    return(rlt)
  }


# Biomarker Selection With WGBS databset and Refine by small-set Plamsa data(5 Paired Tissue-Plasma colon and lung samples and 5 normal plasma)
OneTimeTissuePredictionTrainStage<-function(data,bio){
  input<-data[match(rownames(bio),rownames(data)),]
  Num<-c()
  for(j in seq(0,0.9,0.05)){
    counts1<-apply(input[,grep(".6P|X6.P|X6.T",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts2<-apply(input[,grep(".7P|X7.P|X7.T",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    counts3<-apply(input[,grep("NC.P",colnames(input))],2,function(x) tapply(x,bio[,4],function(x) sum(x>j,na.rm=T)))
    num<-data.frame(id=j,c1=sum(apply(counts1,2,function(x) which.max(x)==2)),
                    c2=sum(apply(counts2,2,function(x) which.max(x)==6)),
                    c3=sum(apply(counts3,2,function(x) which.max(x)==10)))
    Num<-rbind(Num,num)
  }
  print(counts1)
  acc<-sweep(Num, 2,c(1,10,10,5) , `/`)
  rlt<-acc[which.max(rowSums(acc[,2:4])),]
  return(rlt)
}
  
  
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
# feature reduction (WGBS+RRBS missing<60, each plasma category(colon,lung and normal)<60%, Plasma missing<50%)
rm<-frm(data)
data<-data[-rm,]
test=data[,c(grep("X7.T|X6.T",colnames(data)),grep("NC.P",colnames(data))[1:5],grep("X6.P",colnames(data))[c(1,3,4,5,6)],grep("X7.P",colnames(data))[c(1,3,4,5,6)])]
OneTimeTissuePredictionTrainStage(test,bio)
# Biomarker validation in extend cancer and normal plasmas

