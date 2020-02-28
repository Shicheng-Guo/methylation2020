install.packages("KMsurv")

library(survival)
library(KMsurv)
data(tongue); 
attach(tongue);
tongue

my.surv <- Surv(time[type==1], delta[type==1])
fit<-survfit(my.surv~1)
fit
summary(fit)
plot(fit, main="Kaplan-Meier estimate with 95% confidence bounds",xlab="time", ylab="survival function")
