library("fpc")
sampleinfo<-read.table("bam_follow_patient_sample.txt",sep="\t",head=T)
x1<-sampleinfo$Patients
x1
load("fpca.result.RData")  
x2<-array(substr(rownames(fpca.rnaseq.rlt),1,12))    #get the sample order of the cluster analysis
group0<-array(data.matrix(data.frame(sampleinfo[match(x2,x1),6])))  # get all the patients info according to the order of cluster analysis
group1<-group0[-which(is.na(group0))]   # remove the patients who info is "NA"

d <- dist(fpca.rnaseq.rlt, method = "minkowski") # distance matrix
fit <- hclust(d, method="ward") # plot(fit) # display dendogram
group2 <- array(cutree(fit, k=4))[-which(is.na(group0))] # cut tree into 5 clusters

cluster.stats(d,group1 ,group2,compareonly=TRUE) 

fpca.rnaseq.rlt2<-fpca.rnaseq.rlt[-which(is.na(group0)),]
source("heatmap2.R")

gheatmap(fpca.rnaseq.rlt2[,-1],rowname=group1,17)

#dim(fpca.rnaseq.rlt2)
#length(group1)
#kmeans(fpca.rnaseq.rlt2[,-1])
#simu1<-sample(1:3, length(group1), replace = TRUE)
#simu2<-sample(1:4, length(group1), replace = TRUE)
#cluster.stats(d,simu1 ,simu2,compareonly=TRUE) 
#m1<-c(1,1,1,1,2,2,2,2,3,3,3,3)
#m2<-c(1,1,1,1,2,2,2,2,2,2,2,2)
#cluster.stats(d,m1,m2,compareonly=TRUE) 

