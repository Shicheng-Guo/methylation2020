setwd("/home/gsc/Dropbox/Project/APCmeta/R/tcga");
load("LUSC.Rdata");
  type<-array();
  type[which(LUSC4[,13]=="01")]=1;
  type[which(LUSC4[,13]=="11")]=0;
  LUSC5<-cbind(LUSC4,type);
for(i in 1:dim(LUSC5)[1]){for(j in 8:12){
  if((!is.na(LUSC5[i,j])&LUSC5[i,j]<0.3)) LUSC5[i,j]<-0
  else if((!is.na(LUSC5[i,j])&LUSC5[i,j]>0.3)) LUSC5[i,j]<-1
  
}};
  names(LUSC5)[8:12]<-c("x2","x3","x4","x5","x6");
  library(epicalc);
  model_LUSC<-glm(type~x2+x3+x4+x5+x6,data=LUSC5,family=binomial);
  auc_LUSC<-lroc(model_LUSC,title=T)$auc;
 auc_LUSC;

