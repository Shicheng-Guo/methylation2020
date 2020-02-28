
table.binomial.test<-function(region,event.time){
  region1<-table(region)
  region2<-unique(region)
  p=length(region)/(length(region2)*event.time)
  pvalue<-c()
  for(i in 1:length(table(region))){
    ptmp<-binom.test(table(region)[i],event.time,p,alternative="two.sided")$p.value
    pvalue<-c(pvalue,ptmp)
  }
  return(pvalue)
  print(paste("there are",sum(pvalue<0.05/length(pvalue)),"regions were over-preferred in the samples"),sep=" ")
}



setwd("/home/sguo/methylation")
file=list.files(pattern="*mh.cor.RData")
load(file[i])
cancer=substr(file[1],1,4)
edge=0.6
newcor<-cor[which(cor[,2]>edge),]
bed1<-cor2bed(rownames(newcor))
print(nrow(bed1))
number<-c()
region<-c(rownames(newcor))
for(i in 2:length(file)){
  load(file[i])
  cancer=substr(file[i],1,4)
  newcor<-cor[which(cor[,2]>edge),]
  bed2<-cor2bed(rownames(newcor))
  bed1<-Rbedtools(functionstring="intersectBed",bed1,bed2,opt.string="-wa -u")
  ntmp<-c(nrow(bed1),nrow(bed2))
  number<-c(number,nrow(bed1))
  region<-c(region,rownames(newcor))
}
regionx1<-unique(region)
regionx2<-cor2bed(regionx1)

rlt<-table.binomial.test(regionx1)

p<-c()
for(i in 1:11){
p<-c(p,2421*0.4^i)
}

x<-1:11
expectation<-round(p)
expectation<-c(968,387,155,62,25,10,4,2,1,0,0)
observe<-as.numeric(c(602,402,255,215,140,126,138 ,95,198,100,150))
chisq.test(expectation,observe,simulate.p.value = T,B=200000)

chisq.test(matrix(c(602,968,1819,645),2,2,byrow=T))


plot(x,y=expectation,lwd=4,xlim=c(2,11.5),lty=3,type="l",xlab="Number of cancer types",col=2,ylab="Counts")
lines(x,observe,lwd=4,lty=1,col=3)
legend("topright",legend=c("Expectation","Observe"),col=c(2,3),lwd=4,lty=c(3,1),bty = "n")
plot(density(table(table(v))))

