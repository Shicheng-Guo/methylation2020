# This is used to conduct p-value analysis for lung cancer methylation 27k dataset in TCGA(LUAD+LUSC)

setwd("/home/sguo/Downloads/meth27lungcancer")
file<-list.files(pattern="*.txt")
file

data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=T,skip=1,sep="\t",row.names=1)
  data<-cbind(data,tmp[,1])
  print(i)
}
sample<-as.vector(sapply(file,function(x){substr(x,45,56)}))
sample
type<-as.vector(sapply(file,function(x){ substr(x,58,59)}))
colnames(data)<-sample
rownames(data)<-rownames(tmp)
data[1:5,1:5]
save(data,file="/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/data/Meth27KLungCacerTCGA.Rdata")
load("Meth27KLungCacerTCGA.Rdata")

saminfo<-read.table("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/data/saminfo.txt",head=T,sep="\t")
head(saminfo)
samp<-saminfo[match(sample[which(sample %in% saminfo[,1])],saminfo[,1]),]
table(samp$race)

target<-c("cg20530314","cg15343119","cg03567830","cg05722918","cg20881888") # AGTR1, GALR1, NTSR1, SLC5A8, ZMYND10
match<-match(target, rownames(tmp))
data<-data[match,]
data[1:5,1:5]

dat<-data+matrix(rnorm(length(data),0.00001,0.00001),dim(data)[1],dim(data)[2])  # row is variable, col is individuale
data<-dat
# obtain x,y index
output<-matrix(NA,dim(data)[1],5)   # set output matrix ()
x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive
x2<-which(type==names(table(type))[2])   # type 2, normal or resistant
table(type)
for(i in 1:dim(data)[1]){
  a<-try(shapiro.test(as.numeric(data[i,x1]))$p.value)
  b<-try(shapiro.test(as.numeric(data[i,x2]))$p.value)
  if(any(c(a<0.05,b<0.05)==T)){
    tmp1<-try(wilcox.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=F))
    output[i,1]<-tmp1$p.value
    output[i,2]<-mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2]))
    output[i,3]<-"wilcox"
    output[i,4]<-mean(as.numeric(data[i,x1]))
    output[i,5]<-mean(as.numeric(data[i,x2]))
    if(tmp1$p.value<0.05/dim(data)[1]){
      filename=paste(data[i,3],data[i,2],"pdf",sep=".")
      #pdf(filename)
      #type=rep(c("cancer","normal"),49)
      #beeswarm(as.numeric(data[i,4:101])~ type, data = data, method = 'swarm',pch = 16,xlab = '', ylab = 'Methylation Level (Beta)',l
      #boxplot( as.numeric(data[i,4:101])~ type, data = data, add = T,names = c("",""), col="#0000ff22") 
      #dev.off() 
    }
  }else{
    tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2],paired=F)))
    output[i,1]<-tmp1$p.value
    output[i,2]<-tmp1$statistic
    output[i,3]<-"ttest"
    output[i,4]<-mean(as.numeric(data[i,x1]))
    output[i,5]<-mean(as.numeric(data[i,x2]))
    
    if(tmp1$p.value<0.05/dim(data)[1]){
      filename=paste(data[i,3],data[i,2],"pdf",sep=".")
      #pdf(filename)
      #type=rep(c("cancer","normal"),49)
      #beeswarm(as.numeric(data[i,4:101])~ type, data = data, method = 'swarm',pch = 16,xlab = ''
      #boxplot( as.numeric(data[i,4:101])~ type, data = data, add = T,names = c("",""), col="#000
      #dev.off() 
    }
  }
  
  print(i)
}

colnames(output)<-c("Pvalue1","t-statistic","test","Mean1","Mean2")
rownames(output)<-target
write.table(output,file="array.test.txt",sep="\t")
getwd()

type=abs(as.numeric(as.factor(type))-2)
type
dat2<-as.data.frame(cbind(type,t(data)))
head(dat2)
tmp[match(target,rownames(tmp)),2]
data[1:5,1:5]

fit<-glm(type~.,data=dat2,binomial(link = "logit"))
summary(fit) # display results
confint(fit) # 95% CI for the coefficients
exp(coef(fit)) # exponentiated coefficients
exp(confint(fit)) # 95% CI for exponentiated coefficients
prob<-predict(fit, type="response") # predicted values
residuals(fit, type="deviance") # residuals
library(pROC)
g <- roc(type ~ prob, data = dat2)
g
sen<-g$sensitivities[which.max(g$sensitivities+g$specificities)]
spe<-g$specificities[which.max(g$sensitivities+g$specificities)]
sen
spe
plot(g)


auc<-sen<-spe<-c()
for(i in 1:5){
  fit<-glm(type~dat2[,i+1],data=dat2,binomial(link = "logit"))
  prob<-predict(fit, type="response") # predicted values
  g <- roc(type ~ prob, data = dat2)
  sen[i]<-g$sensitivities[which.max(g$sensitivities+g$specificities)]
  spe[i]<-g$specificities[which.max(g$sensitivities+g$specificities)]
  auc[i]<-g$auc # predicted values  
}
colnames(dat2)
sen
spe
auc

write.table(round(cbind(sen,spe,auc),4),file="array.sen.spe.acc.txt",sep="\t",quote=F)


setwd("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/data")
load("Meth27KLungCacerTCGA.Rdata")
cancerresearch<-c("TERT","CDKN2A","WT1","RASSF1")
match<-match(cancerresearch,tmp[,2])
match
dat2<-data[match,]
type=abs(as.numeric(as.factor(type))-2)
dat2<-as.data.frame(cbind(type,t(dat2)))
head(dat2)




setwd("/home/sguo/Dropbox/Project/methylation/Genesky_project/Phase1/result/methyClassification")
data<-read.table("data.txt",head=T,sep="\t",as.is=F,row.names=1)
head(data)
data<-RawNARemove(data)
dat<-data[,c(2,3,4,5,10,11,12,13,18,20)]
head(dat)
#dat<-data[,2:19]  # 1) all the methylation site
type<-abs(as.numeric(as.factor(dat[,10]))-2)
type
rlt<-impute.knn(data.matrix(dat[,1:(dim(dat)[2]-1)]),k = 3)
rlt<-data.frame(rlt$data)

# glasso cross validation
dat2<-cbind(rlt,type)
head(dat2)

auc<-sen<-spe<-c()
col<-c(4,5,7,8,9)
for(i in 1:5){
  i=1
  fit<-glm(type~Gender+Age+Smoking+GALR1+AGTR1+SLC5A8+ZMYND10+NTSR1,data=dat2,binomial(link = "logit"))   # AGTR1"   "GALR1"   "LINE.1"  "NTSR1"   "SLC5A8"  "ZMYND10" "type"  
#  fit<-glm(type~Gender+Age+Smoking+dat2[,col[i]],data=dat2,binomial(link = "logit"))
  prob<-predict(fit, type="response") # predicted values
  g <- roc(type ~ prob, data = dat2)
  sen[i]<-g$sensitivities[which.max(g$sensitivities+g$specificities)]
  spe[i]<-g$specificities[which.max(g$sensitivities+g$specificities)]
  auc[i]<-g$auc # predicted values  
  i
  i=i+1
}
colnames(dat2)
sen
spe
auc
output<-round(cbind(sen,spe,auc),4)
output
write.table(round(cbind(sen,spe,auc),4),file="genesky.sen.spe.acc.txt",sep="\t",quote=F)

pdf("model.expansion.pdf")
plot(sen,type="n",ylim=c(0.4,0.9),ylab="Value",xlab="Biomarker",axes = F,)
axis(1, at = c(1:5), lab = c("GALR1","AGTR1","SLC5A8","ZMYND10","NTSR1"), col = 'transparent', col.tick = 'black', cex.axis = 0.8)
axis(2, at = seq(0.4,0.9,by=0.2), col = 'transparent', col.tick = 'black', cex.axis = 0.8)
lines(sen,type="o",col="red",lty=1)
lines(spe,type="o",col="green",lty=2)
lines(auc,type="o",col="blue",lty=3)
legend("bottomrigh",legend=c("Sensitivity","Specificity","AUC"),lty=c(1:3),col=c("red","green","blue"))
box()
dev.off()


polygon(rep(y[1:2], each = 2), c(0, 28, 28, 0), col = 'lightyellow', border= 'lightyellow')
polygon(c(200, 200, y[1], y[1]), c(0, 28, 28, 0), col = 'lightgray', border= 'lightgray')
matplot(add = T, xlim = c(0, 2667), cbind(rowMeans(xx1),sen$ci1,sen$ci2,rowMeans(xx2),res$ci1,res$ci2),type='l',lty=c(1,4,4,1,4,4),col=c(2,2,2,3,3,3),lwd=c(3,2,2,3,2,2))
legend(x=2000,y=27,box.col="transparent",legend=c("resistant observation","resistant 95% CI","sensitive observation","sensitive 95% CI"),cex=0.8,lty=c(1,4,1,4), bg = "transparent",col=c(2,2,3,3),lwd=c(3,2,3,2))



