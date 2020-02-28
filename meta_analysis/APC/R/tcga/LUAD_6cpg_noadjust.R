setwd("/home/gsc/Dropbox/Project/APCmeta/R/tcga");
load("LUAD.Rdata");
type<-array();
type[which(LUAD4[,13]=="01")]=1;
type[which(LUAD4[,13]=="11")]=0;
LUAD5<-cbind(LUAD4,type);
for(i in 1:dim(LUAD5)[1]){for(j in 8:12){
  if((!is.na(LUAD5[i,j])&LUAD5[i,j]<0.3)) LUAD5[i,j]<-0
  else if((!is.na(LUAD5[i,j])&LUAD5[i,j]>0.3)) LUAD5[i,j]<-1
  
}};
names(LUAD5)[8:12]<-c("x2","x3","x4","x5","x6");
library(epicalc);
model_LUAD<-glm(type~x2+x3+x4+x5+x6,data=LUAD5,family=binomial);
auc_LUAD<-lroc(model_LUAD,title=T)$auc;
auc_LUAD;