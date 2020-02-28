y<-c(rep("A",30),rep("B",75))
x1<-c(abs(rnorm(30,0.7,0.5)),abs(rnorm(75,0.1,0.5)))
x2<-c(abs(rnorm(30,0.1,0.5)),abs(rnorm(75,0.7,0.5)))
x3<-c(abs(rnorm(30,0.2,0.5)),abs(rnorm(75,0.6,0.5)))
t.test(x1[1:30],x1[31:105])
t.test(x2[1:30],x2[31:61])
data=data.frame(y,x1,x2,x3)
ggplot(data,aes(y,x1))+geom_boxplot()+geom_jitter(size=3,colour="red")
ggplot(data,aes(y,x2))+geom_boxplot()+geom_jitter(size=3,colour="blue")
ggplot(data,aes(y,x3))+geom_boxplot()+geom_jitter(size=3,colour="green")

library("ggplot2")
library("randomForest")
t.test(x1[1:30],x1[31:105])
data=data.frame(y,x1,x2)
fit <- randomForest(y~x1+x2,data=data[1:60,])
fit

ggplot(data,aes(x1,x2),col=as.factor(y))+geom_point(aes(colour=factor(y)),size=3)
