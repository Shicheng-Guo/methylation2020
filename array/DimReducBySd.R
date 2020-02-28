#!/usr/bin/R
setwd("/home/gsc/Dropbox/Project/methylation/array")
file<-list.files(pattern="funcPCA_score.txt")
file
for (i in 1:length(file)){
dat<-read.table(file[i],head=T,sep="\t",row.names=1)
sd2<-apply(dat,2, sd)
plot(hist(sd2))
}


