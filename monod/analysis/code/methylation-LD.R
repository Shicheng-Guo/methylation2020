

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
  rlt$r<-as.numeric(r)
  rlt$nobs<-sum(table)
  rlt$phi<-phi
  rlt$chisq<-as.numeric(chisq)
  return(rlt)
}
X<-c()
Y<-c()
for(j in 1:100){
mlc<-c()
Vector<-c()
a<-round(runif(1,1,100))
b<-round(runif(1,1,100))
for(i in 1:100){
  vector<-sample(c(rep("CC",a),rep("CT",b),rep("TC",b),rep("TT",a)),100,replace=T)
  Vector<-c(Vector,vector)
  A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
  A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
  tmp<-data.frame(A1,A2)
  mlc<-rbind(mlc,tmp)
}
x<-LD(Vector)
y<-cor.test(mlc[,1],mlc[,2])
X<-c(X,x$r)
Y<-c(Y,abs(as.numeric(y$estimate)))
print(j)
}
lm(Y~X)
pdf("LDvsCOR.pdf")
plot(x=abs(X),y=Y,xlab="Absolute LD (phi)",ylab="Absolute pearson correlation coefficient (r)",pch=16,col="blue")
abline(a=0,b=1,lwd=3,col="red",lty=6)
abline(h=0.5,lty=4,col="green")
abline(v=0.5,lty=4,col="green")
dev.off()

LD(vector)

data=data.frame(X,Y)
library("ggplot2")
ggplot(data, aes(X, Y)) +  geom_point() +  geom_smooth()



vector<-sample(c(rep("CC",10),rep("CT",100),rep("TC",100),rep("TT",10)),100,replace=T)
table(vector)
LD(vector)
vector<-sample(c(rep("CC",1),rep("CT",1),rep("TC",1),rep("TT",1)),100,replace=T)
LD(vector)
vector<-sample(c(rep("CC",1),rep("CT",2),rep("TC",2),rep("TT",1)),100,replace=T)


? abline
# how to calculate phi
library(vcd)
assocstats(table)
