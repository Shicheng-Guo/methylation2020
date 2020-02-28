
data<-read.table("/home/shg047/oasis/mouse/mergeHapinfo/hapinfo/mf.output.txt",head=T,row.names=1,as.is=T)
phen<-read.table("/home/shg047/oasis/mouse/mergeHapinfo/config.txt")

SRX<-unlist(lapply(colnames(data),function(x) unlist(strsplit(x,"[.]"))[1]))
Type<-as.character(phen[match(SRX,phen[,1]),2])

x1<-which(Type=="mESC")
x2<-which(Type=="Adult")

rlt<-ttestFunction(data,x1,x2)
rlt<-rlt[order(rlt[,3]),]
subset<-subset(rlt,rlt[,3]<0.05 & abs(as.numeric(as.character(rlt[,4])))>0.2)
dim(subset)
write.table(rlt,file="mESC.vs.Adult.5me.differential.txt",sep="\t",quote=F,col.names=T,row.names=F)

ttestFunction<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],5)
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),sum(! is.na(data[i,x1]))<2,sum(! is.na(data[i,x2]))<2,all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]), na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-as.numeric((mean(as.numeric(data[i,x1]),na.rm=T)-mean(as.numeric(data[i,x2]),na.rm=T)))
      output[i,3]<-"t-test"
      output[i,4]<-mean(as.numeric(data[i,x1]),na.rm=T)
      output[i,5]<-mean(as.numeric(data[i,x2]),na.rm=T)
      print(i)
    }
  }
  out<-cbind(rownames(data),output)
  P.Adj<-p.adjust(output[,1],method="fdr")
  out<-data.frame(out[,1],out[,2],P.Adj,out[,3:6])
  colnames(out)=c("RowName","P-value","P.FDR","Delta","Test","M(G1)","M(G2)")
  return(out)
}
