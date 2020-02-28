# Supplementary Figure 11-a
load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")

# add CT specifici MHL
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

# build updated biomarkers
biomark1<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
biomark2<-rbind(biomark1,udpaterlt)

input<-data[match(rownames(biomark2),rownames(data)),]
# automatically select best threshold with 5-fold cross-validation for colon plasma/lung cancer/normal plasma together.
count<-apply(input[,grep("NC.P",colnames(data))],2,function(x) tapply(x,bio$V5,function(x) sum(x>0.02,na.rm=T)))
pdf("supplementaryFigure11.pdf")
par(mfrow=c(4,4))
for(i in 1:nrow(count)){
  hist(count[i,],main=rownames(count)[i],col="darkblue",breaks=30,xlim=c(0,60))
}
dev.off()

# Supplementary Figure 11-b
par(mfrow=c(2,2))
hist(CCP[2,],breaks=30,xlim=c(0,60),col="darkblue",border="darkblue",xlab="Colon",ylab="Counts of cs-MHL",main="Colon")
hist(CCP[11,],breaks=10,xlim=c(0,60),col="darkblue",border="darkblue",xlab="CT",ylab="Counts of cs-MHL",main="Cancer Tissue")
dim(normal)



