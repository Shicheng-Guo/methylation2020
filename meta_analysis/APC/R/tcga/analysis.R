setwd("/home/gsc/Dropbox/Project/APCmeta/R/tcga")
load("LUAD.Rdata")
load("LUSC.Rdata")
luad<-LUAD4
lusc<-LUSC4




#for correlation of each cpg site
library("corrplot")
cor1<-matrix(NA,6,6)
for (i in 1:6){
  for (j in 1:6){
  cor1[i,j]<-cor(luad[,i+6],luad[,j+6],use="complete.obs")
  }
}

cor2<-matrix(NA,6,6)
for (i in 1:6){
  for (j in 1:6){
    cor2[i,j]<-cor(lusc[,i+6],lusc[,j+6],use="complete.obs")
  }
}

pdf("corrplot.pdf")
par<-par(mfrow=c(1,2))
colnames(cor1)<-colnames(luad)[7:12]
colnames(cor2)<-colnames(luad)[7:12]
rownames(cor1)<-colnames(luad)[7:12]
rownames(cor2)<-colnames(luad)[7:12]
corrplot(cor1,method="circle")
corrplot(cor2,method="circle")
dev.off()

write.table(cor1,file="cor1.txt",sep="\t",col.names=colnames(luad)[7:12],row.names=colnames(luad)[7:12])
write.table(cor2,file="cor2.txt",sep="\t",col.names=colnames(luad)[7:12],row.names=colnames(luad)[7:12])



dim(lusc[which(lusc$casecontrol=="01"),])
dim(luad[which(luad$casecontrol=="01"),])


# to identify normal stage distribution
x1<-(luad[which(luad$casecontrol=="01"),])
x2<-(luad[which(luad$casecontrol=="11"),])
x3<-(lusc[which(lusc$casecontrol=="01"),])
x4<-(lusc[which(lusc$casecontrol=="11"),])


t1<-table(x1$stage)
t2<-table(x2$stage)
t3<-table(x3$stage)
t4<-table(x4$stage)
m<-rbind(t1,t2,t3,t4)
write.table(m,file="stage.txt",sep="\t")
getwd()

s1<-lusc$age[which(lusc$casecontrol=="01")]
s2<-lusc$age[which(lusc$casecontrol=="11")]
s3<-luad$age[which(luad$casecontrol=="01")]
s4<-luad$age[which(luad$casecontrol=="11")]

length(which(s1>65))
length(which(s1<=65))
range(s1,na.rm=T)

length(which(s3>65))
length(which(s3<=65))
range(s3,na.rm=T)


paste(round(mean(s1,na.rm=T),2)," (",round(sd(s1,na.rm=T),2),")",sep="")
paste(round(mean(s2,na.rm=T),2)," (",round(sd(s2,na.rm=T),2),")",sep="")
paste(round(mean(s3,na.rm=T),2)," (",round(sd(s3,na.rm=T),2),")",sep="")
paste(round(mean(s4,na.rm=T),2)," (",round(sd(s4,na.rm=T),2),")",sep="")


s1<-table(lusc$gender[which(lusc$casecontrol=="01")])
s2<-table(lusc$gender[which(lusc$casecontrol=="11")])

s3<-table(luad$gender[which(luad$casecontrol=="01")])
s4<-table(luad$gender[which(luad$casecontrol=="11")])

sum(luad$casecontrol=="11")

M<-as.table(rbind(s1,s2))
N<-as.table(rbind(s3,s4))


s1<-table(lusc$stage[which(lusc$casecontrol=="01")])
s3<-table(luad$stage[which(luad$casecontrol=="01")])


#There is no significant difference of the age and gender between case and control (S.table). 
library(epicalc)


luad$casecontrol<-as.numeric(as.character(luad$casecontrol))
luad$casecontrol[which(luad$casecontrol==11)]<-0
lusc$casecontrol<-as.numeric(as.character(lusc$casecontrol))
lusc$casecontrol[which(lusc$casecontrol==11)]<-0


model1<-step(glm(casecontrol~cg01240931+cg15020645+cg16970232+cg21634602+cg20311501+cg24332422,data=luad,family=binomial))
predict(model1,type="response")
lroc(model1, graph = TRUE, add = FALSE, title = FALSE,  auc.coords = NULL)
summary(model1)

model2<-glm(casecontrol~cg01240931+cg15020645+cg16970232+cg21634602+cg20311501+cg24332422,data=lusc,family=binomial)
lroc(model2, graph = TRUE, add = FALSE, title = FALSE,  auc.coords = NULL)
predict(model2,type="response")
summary(model2)

roc.area(model)
lroc(model, graph = TRUE, add = FALSE, title = FALSE,  auc.coords = NULL)


cgsite<-c("cg01240931","cg15020645","cg16970232","cg21634602","cg20311501","cg24332422")
for (i in 1:6){
  model1<-(glm(casecontrol~cg01240931,data=luad,family=binomial))
  summary(model1)
  model2<-(glm(casecontrol~cg01240931,data=lusc,family=binomial))
  summary(model2)
}

model2<-glm(casecontrol~cg01240931+age+gender,data=lusc,family=binomial)
summary(model2)



names(model.frame(model))
model.frame(model)[,1][,2:1]

names(model)

head(table(model$fitted.value,model$y))




