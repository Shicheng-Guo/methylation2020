
mhl<-c()
for(i in seq(5,1000,by=1)){
w=4:i
p=runif(length(w),0,0.3)
mhl1=sum(p*w)/sum(w)
mhl<-c(mhl,mhl1)
}
plot(mhl,ylab="MHL",xlab="Number of Methylation Haplotypes",col="red",type="l",lwd=3)

qqplotR<-function(p){
  observed <- sort(p)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed)) 
  lexp <- -(log10(expected / (length(expected)+1)))
  data.frame(lobs,lexp)
}

p1<-c()
z<-rnorm(1000,0,1)
for(i in 1:1000000){
x<-rnorm(1000,0,1)
y<-rnorm(1000,0,1)
fit<-lm(z~x+y+x*y)
tmp<-summary(fit)$coefficients[4,4]
p1<-c(tmp,p1)
}
p2<-c()
r<-rnorm(10000,0,0.1)
for(i in 1:10000000){
  x<-rnorm(1000,0,1)
  y<-rnorm(1000,0,1)
  z<-x+y+r*x*y
  fit<-lm(z~x+y+x*y)
  tmp<-summary(fit)$coefficients[4,4]
  p2<-c(tmp,p2)
}

p1<-sample(p1,10000)
p2<-sample(p2,10000)

plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
pt<-seq(0,0.2,by=0.05)
for(i in 1:length(pt)){
a<-pt[i]
p<-c(sample(p1,round((1-a)*length(p1))),sample(p2,round(a*length(p1))))
pp1<-qqplotR(p)
points(pp1$lexp,pp1$lobs, pch=1, cex=.4, col=rainbow(7)[i+1])   # PC3
legend("bottomright",legend=paste("Postive%=",pt,sep=""),pch=1,col=c(rainbow(7)[2:(length(pt)+1)]),bty="n")
}








