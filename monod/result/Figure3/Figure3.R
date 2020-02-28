setwd("/home/sguo/Dropbox/Project/methylation/monod/result/Figure3")
rm(list=ls())
load("wgbs.mhl.heatmap.RData")
dim(mydata)

gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,mean)))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(data[i,][which(index==group[j])]))/10^(mean(data[i,][which(index==gmax)])))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    gsi<-c(gsi,gsit)
    #  print(c(gmax,gsit))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi)
  return(rlt)
}
TopGSIByCategory<-function(gsi,top=150){
  GSIRlt<-c()
  group<-names(table(gsi$group))
  rank<-c(rep(top,length(group)))
  for (i in 1:length(group)){
    subset=gsi[which(gsi$group==group[i]),]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
  }
  return(na.omit(GSIRlt))
}
gsirlt<-gsi(mydata)
topgsi<-TopGSIByCategory(gsirlt)
input<-mydata[match(topgsi[,1],rownames(mydata)),]
mydata <- scale(input) # standardize variables
d <- dist(t(mydata), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
fit$labels[fit$order]




