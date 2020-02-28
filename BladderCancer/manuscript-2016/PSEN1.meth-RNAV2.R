
# Analysis the correlation of methylation and expression for PSEN1
logit<-function (p, percents = range.p[2] > 1, adjust) {
  range.p <- range(p, na.rm = TRUE)
  if (percents) {
    if (range.p[1] < 0 || range.p[1] > 100) 
      stop("p must be in the range 0 to 100")
    p <- p/100
    range.p <- range.p/100
  }
  else if (range.p[1] < 0 || range.p[1] > 1) 
    stop("p must be in the range 0 to 1")
  a <- if (missing(adjust)) {
    if (isTRUE(all.equal(range.p[1], 0)) || isTRUE(all.equal(range.p[2], 
                                                             1))) 
      0.025
    else 0
  }
  else adjust
  if (missing(adjust) && a != 0) 
    warning(paste("proportions remapped to (", a, ", ", 1 - 
                    a, ")", sep = ""))
  a <- 1 - 2 * a
  log((0.5 + a * (p - 0.5))/(1 - (0.5 + a * (p - 0.5))))
}
input=("C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\bladdercancer\\manuscript-2016\\PSEN1.Meth.RNAV2.txt")
data=read.table(input,sep="\t")
plot(data[,2],data[,4])
input=data[,c(1,2,4)]
colnames(input)=c("ID","Meth","EXP")

# linear model
input$Meth=logit(input$Meth)
input$EXP=log(input$EXP,2)
plot(density(na.omit(input$Meth)))

input$EXP=log(input$EXP,2)
x<-na.omit(input$Meth)
x=logit(x)
x<-(x-mean(x))/sd(x)
input$Meth<-x
input<-na.omit(input)
cor(input$EXP,input$Meth)
? cor
plot(input$EXP~input$Meth,xlab="Z-score ( logit( Beta ) )", ylab="log2(Expression)",col="blue",main="PSEN1",cex=1.5,pch=16)

fit<-lm(input$EXP~input$Meth)
abline(fit,lwd=2,col="red")

summary(fit)
plot(fit)
# linear mixture model
library(nlme)
input=na.omit(input)
m1.nlme = lme(EXP~ Meth,random = ~ 1|ID,data = input)
summary(m1.nlme)
plot(ranef(m1.nlme))
plot((m1.nlme))

anova(m1.nlme)

fm1 <- lme(distance ~ age, data = Orthodont) # random is ~ age
fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
summary(fm1)
summary(fm2)
plot(ranef(fm2))
plot((fm2))


