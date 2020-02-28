setwd("/home/gsc/Dropbox/Project/methylation/Genesky_project/Phase1/result")
data<-read.table("data.txt",head=T,sep="\t",as.is=F)
pheno=read.table("pheno.txt",head=T,sep="\t",as.is=F)
dat<-merge(pheno,data,by="IID")
phe<-abs(as.numeric(dat$phen)-2)


# phenotype distribution



# Pvalue based on Logistic regression

out<-data.frame()
for(i in 10:dim(dat)[2]){
sum.coef<-summary(glm(phe~dat[,i]+Gender+Age+Smoking,data=dat,family=binomial(link = "logit"),control = list(maxit = 50)))
lower.ci<-log(exp(sum.coef$coefficients[,1]+1.96*sum.coef$coefficients[,2]),10)
upper.ci<-log(exp(sum.coef$coefficients[,1]-1.96*sum.coef$coefficients[,2]),10)
fit<-round(log(exp(sum.coef$coefficients[,1])[2],base=10),2)
low<-round(lower.ci[2],2)
up<-round(upper.ci[2],2)
pvalue<-format(sum.coef$coefficients[2,4],scientific=T, digits=3)
out[i-9,1]<-paste(fit," (",up,"-",low,")",sep="")
out[i-9,2]<-pvalue
}
out<-matrix(unlist(out),15,2)
write.table(out, file="glm.result.txt",quote=F,sep="\t")

# Pvlaue based on paired t-test
data<-read.table("data.txt",head=T,sep="\t",as.is=F,row.names=1)
ind<-dim(data)[1]
type=rep(c("Case","Control"),150)
t<-try(lapply(data,function(x) t.test(x[seq(1,ind,by=2)],x[seq(2,ind,by=2)],paired=T)))
pvalue<-lapply(t, function(x) x$p.value)
tvalue<-lapply(t, function(x) x$statistic)
xmean<-lapply(t, function(x) x$estimate[1])
ymean<-lapply(t, function(x) x$estimate[2])

mean<-cbind(xmean,ymean)
write.table(mean,file="mean.txt",sep="\t",col.names=NA, row.names=T)

pvalueadj<-p.adjust(pvalue, method = p.adjust.methods, n = length(pvalue))
res<-rbind(pvalueadj,tvalue)
write.table(res,file="dataresult.txt",sep="\t",col.names=NA, row.names=T)

# beeswarm plot for paired sample
#install.packages("beeswarm")
require(beeswarm)
op<-par(mfrow(3,3))
pdf("beeswarm.pdf")
for (i in 1:15){
beeswarm(data[,i]~ type, data = data, method = 'swarm',pch = 16,xlab = '', ylab = 'Methylation Level',labels = c('Cancer', 'Control'),main=colnames(data)[i])
boxplot( data[,i]~ type,data = data, add = T,names = c("",""), col="#0000ff22") 
}
dev.off()
par(op)



#Table 1. Characteritic of patients
case<-dat[which(dat$pheno=="cancer"),]
con<-dat[which(dat$pheno=="normal"),]
colnames(case)
density(case$Age)
table(case$Gender)
table(case$type)
table(case$TNM)
table(case$differention)

#Table 2. Methylation level of the candidate genes in NSCLC and controls
library("gplots")
colnames(dat)[16]<-"LINE-1"
colnames(dat)[24]<-"Reference"
save(dat,file="dat.RData")
pdf("heatmap.2.pdf")
heatmap.2(cor(dat[,10:24], use="complete.obs"), col = topo.colors(75),key=T,keysize=1, trace="none")
dev.off()

head(dat)
#age associated methylation
rlt_age<-(lapply(dat[10:24],function(x) summary(lm(x~dat$Age))$coefficients[2,4]))
rlt_gender<-(lapply(dat[10:24],function(x) summary(lm(x~dat$Gender))$coefficients[2,4]))
rlt_Smoking<-(lapply(dat[10:24],function(x) summary(lm(x~dat$Smoking))$coefficients[2,4]))
rlt_type<-(lapply(dat[10:24],function(x) summary(lm(x~dat$type))$coefficients[2:3,4]))
rlt_diff<-(lapply(dat[10:24],function(x) summary(lm(x~dat$differention))$coefficients[2,4]))
rlt_tnm<-(lapply(dat[10:24],function(x) summary(lm(x~dat$TNM))$coefficients[2,4]))
# combination
rlt<-(lapply(dat[10:24],function(x) summary(lm(x~dat$Age+dat$Gender))$coefficients[2:3,4]))


cor(as.numeric(dat[,3]), dat[,16],use = "complete.obs")

x<-which(dat$pheno=="normal")
x<-which(dat$pheno=="cancer")

cor(as.numeric(dat[x,3]), dat[x,16],use = "complete.obs")





