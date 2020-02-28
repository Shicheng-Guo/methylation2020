x<-na.omit(x)
y<-na.omit(y)

rlt1<-c()
rlt2<-c()
for(i in c(1,1.5,2,4,9,99,999)){
  Mean<-c()
  Meanzz<-c()
  for(j in 1:50){
    means<-mean(c(sample(x,1000,replace=T),sample(y,1000*i,replace=T)),na.rm=T)
    Mean<-c(Mean,means)
    meanszz<-means+sign(0.5-means)*rnorm(1,0.05,0.005)
    Meanzz<-c(Meanzz,meanszz)
  }
  rlt1<-c(rlt1,mean(Mean)/sd(Mean))
  rlt2<-c(rlt2,mean(Meanzz)/sd(Meanzz))
}
rlt1
rlt2
names(rlt1)=c("50%","40%","30%","20%","10%","1.0%","0.1%")
plot(rlt1,ylim=c(10,150),type="n",ylab="Mean/SD")
lines(rlt1,col="red",lwd=3)
lines(rlt2,col="blue",lwd=3)
