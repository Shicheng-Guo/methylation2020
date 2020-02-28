

file=list.files(pattern="*pairwiseR2$")
for(i in 1:length(file)){
  data<-read.table(file[i],sep="\t",skip=1)
  id<-unique(data[,1])
  gap1<-c()
  gap2<-c()
  r2<-c()
  for(j in 1:length(id)){
    tmp<-data[data[,1]==id[j],]
    if(nrow(tmp)>2 && sd(tmp[,3],na.rm=T) != 0  && summary(lm(tmp[,3]~tmp[,2],na.action="na.omit"))$coefficients[2,1]<0){
      r2<-c(r2,tmp[,3])
      gap1<-c(gap1,tmp[,2])
      gap2<-c(gap2,tmp[,2]/max(tmp[,2]))
      print(j)
      if (length(gap1)>5000) break
    }  
  }
  par(mfrow=c(2,2))
  pdf(paste(file[i],".paisewise.R2.pdf",sep=""))
  plot(r2~gap1)
  plot(r2~gap2)
  dev.off()
}

library(RColorBrewer)
require(KernSmooth)
my.cols <- rev(brewer.pal(11, "RdYlBu"))
file=list.files(pattern="*pairwiseR2$")
for(i in 1:length(file)){
  data<-read.table(file[i],sep="\t",skip=1)
  id<-unique(data[,1])
  gap1<-c()
  gap2<-c()
  r2<-c()
  for(j in 1:length(id)){
    tmp<-data[data[,1]==id[j],]
    if(nrow(tmp)>2 && sd(tmp[,3],na.rm=T) != 0 && summary(lm(tmp[,3]~tmp[,2],na.action="na.omit"))$coefficients[2,1]<0){
      r2<-c(r2,tmp[,3])
      gap1<-c(gap1,tmp[,2])
      gap2<-c(gap2,tmp[,2]/max(tmp[,2]))
      print(j)
      if (length(gap1)>5000) break
    }  
  }
  par(mfrow=c(2,2))
  pdf(paste(file[i],".paisewise.R2.smooth.pdf",sep=""))
  smoothScatter(x=gap1,y=r2, nrpoints=.3*length(gap1), colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1")
  smoothScatter(x=gap2,y=r2, nrpoints=.3*length(gap2), colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1")
  dev.off()
}

#!/home/shg047/software/R-patched/bin/R
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
PairwiseR2plot(R2file=args[1],samplingsize=5000)

PairwiseR2plot <- function(R2file=args[1],samplingsize=5000){
  file=R2file
  library(RColorBrewer)
  require(KernSmooth)
  
  my.cols <- rev(brewer.pal(11, "RdYlBu"))
    data<-read.table(file,sep="\t",skip=1)
    id<-unique(data[,1])
    gap1<-c()
    gap2<-c()
    r2<-c()
  
    for(j in 1:length(id)){
      tmp<-data[data[,1]==id[j],]
      tmp<-na.omit(tmp)
      if(nrow(tmp)>2 && sd(tmp[,3],na.rm=T) != 0 && summary(lm(tmp[,3]~tmp[,2],na.action="na.omit"))$coefficients[2,1]<0){
        r2<-c(r2,tmp[,3])
        gap1<-c(gap1,tmp[,2])
        gap2<-c(gap2,tmp[,2]/max(tmp[,2],na.rm=T))
        print(j)
        if (length(gap1)>samplingsize) break
        }        
  }
  pdf(paste(file,".paisewise.R2.smooth.pdf",sep=""))
  par(mfrow=c(2,2))
  smoothScatter(x=gap1,y=r2, nrpoints=.3*length(gap1), colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1")
  smoothScatter(x=gap2,y=r2, nrpoints=.3*length(gap2), colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1")
  plot(r2~gap1)
  plot(r2~gap2)
  dev.off()  
}



for i in `ls N37*txt.pairwiseR2`
do
R CMD BATCH "--args $i" pairwise.R &
done




