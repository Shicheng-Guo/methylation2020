
LD<-function(vector){
  rlt<-list()
  table<-matrix(table(vector),2,2)
  pAB=table[1,1]/sum(table)
  pA=(2*table[1,1]+table[2,1]+table[1,2])/(2*sum(table))
  pB=(2*table[2,2]+table[2,1]+table[1,2])/(2*sum(table))
  pa=1-pA
  pb=1-pB
  D=pAB-pA*pB
  if(D>0){
    Dmax=min(pA*pb,pa*pB)
  } else{
    Dmax=max(-pA*pB,-pa*pb)
  } 
  Dp=D/Dmax
  r=Dp/sqrt(pA*pa*pB*pb)
  test<-chisq.test(table)
  chisq<-test$statistic
  
  A1<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))
  A2<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))  
  fit<-cor.test(A1,A2)
  corr=fit$estimate
  rlt$nobs<-sum(table)
  rlt$corr=as.numeric(corr)
  rlt$corr.p=fit$p.value
  rlt$Dp=Dp
  rlt$r=r
  rlt$p=test$p.value
  return(rlt)
}


tmpp<-c()
for(j in 1:1000){
  a=abs(rnorm(1,3,10))
  b=abs(rnorm(1,3,10))
  c=abs(rnorm(1,3,10))
  d=abs(rnorm(1,3,10))
  mlc<-c()
  r<-c()
  for(i in 1:100){
    vector<-sample(c(rep("CC",a),rep("CT",b),rep("TC",c),rep("TT",d)),100,replace=T)
    A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
    A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
    tmp<-data.frame(A1,A2)
    mlc<-rbind(mlc,tmp)
    r<-rbind(r,LD(vector))
  }
  
  rpadj<-mean(na.omit(unlist(data.frame(r)[,6])))
  corpadj<-cor.test(mlc[,1],mlc[,2])$p.value
  tmp<-c(rpadj,corpadj)
  tmp
  tmpp<-rbind(tmpp,tmp)
}

######
library("ggplot2")
tmpp<-c()
for(j in 1:10){
  a=abs(rnorm(1,5,10))+2
  b=abs(rnorm(1,5,10))+2
  c=abs(rnorm(1,5,10))+2
  d=abs(rnorm(1,5,10))+2
  mlc<-c()
  r<-c()
  for(i in 1:100){
    vector<-sample(c(rep("CC",a),rep("CT",b),rep("TC",c),rep("TT",d)),100,replace=T)
    A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
    A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
    tmp<-data.frame(A1,A2)
    mlc<-rbind(mlc,tmp)
    r<-rbind(r,LD(vector))
  }
  
  rpadj<-mean(na.omit(unlist(data.frame(r)[,6])))
  corpadj<-cor.test(mlc[,1],mlc[,2])$p.value
  tmp<-c(rpadj,corpadj,sd(na.omit(mlc[,1])))
  tmp
  tmpp<-rbind(tmpp,tmp)
}

colnames(tmpp)<-c("LD","Cor","sd")
tmpp<-data.frame(tmpp)

library("ggplot2")
setwd("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis")
load("LDCor.RData")
f<-ggplot(tmpp,aes(-log(LD,10),-log(Cor,10)))+geom_point(size=6)+ geom_smooth(method = "lm", se = FALSE)
f<-f+xlab("LD:-log(P,10)")+ylab("Cor:-log(P,10)")+theme_bw()+theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
pdf("example.pdf",width=5,height=5)
print(f)
dev.off()

x<--log(tmpp[,1],10)
y<--log(tmpp[,2],10)
xx<-data.frame(x,y)


