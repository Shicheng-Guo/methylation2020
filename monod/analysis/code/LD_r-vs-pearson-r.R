vector<-Vector

LD<-function(vector){
  rlt<-list()
  table<-matrix(table(vector),2,2)
  pAB=table[1,1]/sum(table)
  pA=(table[1,1]+table[1,2])/(sum(table))
  pB=(table[1,1]+table[2,1])/(sum(table))
  pa=1-pA
  pb=1-pB
  D=pAB-pA*pB
  if(D>0){
    Dmax=min(pA*pb,pa*pB)
  } else{
    Dmax=max(-pA*pB,-pa*pb)
  } 
  Dp=D/Dmax
  r=D/sqrt(pA*pa*pB*pb)
  test<-chisq.test(table,correct = F)
  chisq<-test$statistic
  phi=as.numeric(sqrt(test$statistic/length(vector)))
  A1<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))
  A2<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))  
  fit<-cor.test(A1,A2)
  rlt$corr=as.numeric(fit$estimate)
  rlt$corr.p=as.numeric(fit$p.value)
  rlt$p=test$p.value
  rlt$Dp=Dp
  rlt$D=D
  rlt$r<-as.numeric(r)
  rlt$nobs<-sum(table)
  rlt$phi<-phi
  rlt$chisq<-as.numeric(chisq)
  return(rlt)
}
X<-c()
Y<-c()
for(j in 1:500){
  mlc<-c()
  Vector<-c()
  a<-round(runif(1,0,5))
  b<-round(runif(1,0,105))
  for(i in 1:100){
    vector<-sample(c(rep("CC",b),rep("TT",b),rep("CC",a),rep("TT",a)),10,replace=T)
    Vector<-c(Vector,vector)
    A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
    A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
    tmp<-data.frame(A1,A2)
    mlc<-rbind(mlc,tmp)
  }
  x<-LD(Vector)
  table(Vector)
  x
  y<-cor.test(mlc[,1],mlc[,2])
  x
  X<-c(X,x$r)
  Y<-c(Y,abs(as.numeric(y$estimate)))
  print(j)
}

for(j in 1:300){
  mlc<-c()
  Vector<-c()
  ldr<-c()
  a<-round(runif(1,1,100))
  b<-round(runif(1,1,100))
  c<-round(runif(1,1,100))
  d<-round(runif(1,1,100))
  for(i in 1:100){
    vector<-sample(c(rep("CC",a),rep("CT",b),rep("TC",c),rep("TT",d)),100,replace=T)
    Vector<-c(Vector,vector)
    A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
    A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
    tmp<-data.frame(A1,A2)
    mlc<-rbind(mlc,tmp)
    ldr<-c(ldr,LD(vector)$r)
  }
  y<-cor.test(mlc[,1],mlc[,2])
  X<-c(X,median(ldr))
  Y<-c(Y,abs(as.numeric(y$estimate)))
  print(j)
}

pdf("LDvsCOR2.pdf",width=5,height=5)
plot(x=abs(X)^2,y=Y^2,xlim=c(0,1),xlab="LD (r2)",ylab="Pearson correlation coefficient (r2)",pch=16,col="blue")
abline(a=0,b=1,lwd=3,col="red",lty=6)
abline(a=0.25,b=1,lwd=3,col="red",lty=6)
abline(a=-0.25,b=1,lwd=3,col="red",lty=6)
dev.off()




# Answer Dnih's question about imprinting regions

X<-c()
Y<-c()
for(j in 1:100){
  mlc<-c()
  Vector<-c()
  a<-round(runif(1,100,200))
  b<-round(runif(1,1,10))
  for(i in 1:100){
    vector<-sample(c(rep("CC",a),rep("CT",b),rep("TC",b),rep("TT",a)),10,replace=T)
    Vector<-c(Vector,vector)
    A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
    A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
    tmp<-data.frame(A1,A2)
    mlc<-rbind(mlc,tmp)
  }
  x<-LD(Vector)
  y<-cor.test(mlc[,1],mlc[,2])
  x
  X<-c(X,x$r)
  Y<-c(Y,abs(as.numeric(y$estimate)))
  print(j)
}
plot(x=abs(X),y=Y,xlab="Absolute LD (r)",ylab="Absolute pearson correlation coefficient (r)",pch=16,col="blue")
abline(a=0,b=1,lwd=3,col="red",lty=6)
abline(a=0.25,b=1,lwd=3,col="red",lty=6)
abline(a=-0.25,b=1,lwd=3,col="red",lty=6)


