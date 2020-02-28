n1=100;
n2=100;
d1=200;
d2=200;
dat1<-matrix(rnorm(n1*d1,1,1),n1,d1)
dat2<-matrix(rnorm(n2*d2,4,1),n2,d2)
dat = rbind(dat1,dat2)
y = as.factor(c(rep(1,n1),rep(2,n2)))
dwdfit = kdwd(dat,y=y,scaled=T)
w = dwdfit@w[[1]]
dat3 = t(t(dat1) - mean(dat1%*%w)*w)
dat4 = t(t(dat2) - mean(dat2%*%w)*w)
 

