
###   Figure LD plot, heatmap for LD blocks
library("grDevices")
col=colorRampPalette(c("white", "red"))(20)
M <- matrix(runif(400),20,20)
M[lower.tri(M)] <- NA
M[M>0.2]<-1
image(M,col = col,frame=F,xaxt="n",yaxt="n")


###   BSP circle plot base on matrix
M<-matrix(sample(c(0,1,1,1),1000,replace=T),10,10)
col=colorRampPalette(c("white", "red"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)),ylim=c(0,nrow(M)))
for(i in 1:nrow(M)){
  for(j in 1:ncol(M)){
    points(i,j,col=1,pch=circle[M[i,j]+1],cex=1.5)
    
  }
}

###   BSP circle plot base on matrix
M<-matrix(sample(c(1),25,replace=T),5,5)
M[,c(1:2)]<-0
circlePlot(M)

circlePlot<-function(M){
color=c(1,2)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)+1),ylim=c(0,nrow(M)+1),bty="n",xaxt="n",yaxt="n")
for(i in 1:nrow(M)){
  for(j in 1:ncol(M)){
    points(j,i,col=color[M[i,j]+1],pch=19,cex=1.5)
  }
}
}

M<-matrix(sample(c(0,1),25,replace=T),5,5)
color=c(1,2)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)+1),ylim=c(0,nrow(M)+1),bty="n",xaxt="n",yaxt="n")
for(i in 1:nrow(M)){
  for(j in 1:ncol(M)){
    points(j,i,col=color[M[i,j]+1],pch=19,cex=1.5)
  }
}




library("genetics")


g1 <- genotype( c('T/A', NA, 'T/T', NA, 'T/A', NA, 'T/T', 'T/A',
                  'T/T', 'T/T', 'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T',
                  NA, 'T/A', 'T/A', NA) )
g2 <- genotype( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
                  'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
                  'C/A', 'C/A', 'C/A', 'A/A') )
g3 <- genotype( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
                  'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A', 'T/T',
                  'T/A', 'T/A', 'T/A', 'T/T') )
# Compute LD on a single pair
LD(g1,g2)

x<-seq(0,1,0.1)
y<-seq(0,1,0.1)
cor(x,y)
data.frame(x,y)



vector<-sample(c("CC","TT","CT","TC"),100,replace=T)
table<-matrix(table(vector),2,2)
pAB=table[1,1]/sum(table)
pA=(2*table[1,1]+table[2,1]+table[1,2])/(2*sum(table))
pB=(2*table[2,2]+table[2,1]+table[1,2])/(2*sum(table))
D=pAB-pA*pB

LD<-function(table){
  pAB=table[1,1]/sum(table)
  pA=(2*table[1,1]+table[2,1]+table[1,2])/(2*sum(table))
  pB=(2*table[2,2]+table[2,1]+table[1,2])/(2*sum(table))
  pa=1-pA
  pb=1-pB
  D=pAB-pA*pB
  if(D>0){
    Dmax=min(pA*pb,pa*pB)
  } else{
    Dmax=max(???pA*pB,???pa*pb)
  } 
  Dp=D/Dmax
  r=???D/sqrt(pA*pa*pB*pb)
  return(c(Dp,r))
}

LD(table)




