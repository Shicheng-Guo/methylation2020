#######################################################################################################################
###   Title : Methylation haplotype load, absolute methylation level and methylation entropy
###   Lab:  Kun's Lab (http://zhang.openwetware.org/)   
###   Time :  Sep/29/2015 
###   Update: re-define MHL to identify continous long methylation and non-methylation reads
#######################################################################################################################
library("entropy")
methentropy<-function(hapinfo){
  entropy<-entropy(table(hapinfo),unit="log2")/nchar(hapinfo[1])
  return(entropy)
}

hap2methfreq<-function(hapinfo){
  hapinfo<-paste(hapinfo,collapse="")
  hapinfo<-strsplit(hapinfo,"")
  num_C<-length(grep("C",unlist(hapinfo)))
  num_T<-length(grep("T",unlist(hapinfo)))
  methfreq<-  num_C/(num_C+num_T)
  return(methfreq)
}

hap_comb_tmp<-function(hapinfo.single){
  # hap_com.single is a character of haplotype, such as "TCTCTT"
  hap_com.single<-c()
  hapinfo.single<-as.character(hapinfo.single)
  for(i in 1:nchar(hapinfo.single)){
    for(j in i:(nchar(hapinfo.single))){
      hap_com.single<-c(hap_com.single,(substr(hapinfo.single,i,j)))
    } 
  }
  return(hap_com.single)
}

hap_comb<-function(hapinfo){
  # haplotype is array of observed methylation haplotype
  hapinfo<-as.list(hapinfo)
  hap_comb<-unlist(lapply(hapinfo,hap_comb_tmp))
  return(hap_comb)
}

epipoly<-function(x){
  
}

mhl<-function(hap_comb_array){
  # hap_comb_array is combinarial methylation haplotype
  mhl_value<-c()
  hap_comb_array<-unlist(hap_comb_array)
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  unmeth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="T"
  hap_comb_array.CG<-c(hap_comb_array[meth_hap],hap_comb_array[unmeth_hap])
  for(i in unique(nchar(hap_comb_array.CG))){
    input_tmp<-which(nchar(hap_comb_array.CG)==i)  
    nmeth<-length(grep(paste(rep("C",i),collapse=""),hap_comb_array.CG[input_tmp]))
    nunmeth<-length(grep("T",hap_comb_array.CG[input_tmp]))
    mhl_value<-c(mhl_value,i*(nmeth)/sum(nchar(hap_comb_array)==i))
    print(sum(nchar(hap_comb_array)==i))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}


n=100
M<-matrix(sample(c(0,0,0,1),4*n,replace=T),n,4)
M
col=colorRampPalette(c("white", "red"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)),ylim=c(0,nrow(M)))
for(i in 1:nrow(M)){
  for(j in 1:ncol(M)){
    points(i,j,col=1,pch=circle[M[i,j]+1],cex=1.5)
  }
}

dim(M)



###   BSP circle plot base on matrix
## Figure sub1
N=20
C=5
M<-matrix(sample(c(sample(0,1,replace=T),c(sample(1,50,replace=T))),N*C,replace=T),N,C)
col=colorRampPalette(c("white", "red"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)+1),ylim=c(0,nrow(M)+1))
for(i in 1:ncol(M)){
  for(j in 1:nrow(M)){
    points(i,j,col=1,pch=circle[M[j,i]+1],cex=1.5)
  }
}

# plot horizontal line
for(j in 1:nrow(M)){    
  abline(h=j)    
}

## Figure sub2
N=20
C=5
M<-matrix(sample(c(sample(0,1,replace=T),c(sample(1,50,replace=T))),N*C,replace=T),N,C)
col=colorRampPalette(c("white", "red"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)+1),ylim=c(0,nrow(M)+1))
for(i in 1:ncol(M)){
  for(j in 1:nrow(M)){
    points(i,j,col=1,pch=circle[M[j,i]+1],cex=1.5)
  }
}

# plot horizontal line
for(j in 1:nrow(M)){    
    abline(h=j)    
}


## Fgiure sub3


N=100
C=5
M<-matrix(sample(c(sample(0,100,replace=T),c(sample(1,replace=T))),N*C,replace=T),N,C)
col=colorRampPalette(c("white", "red"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)+1),ylim=c(0,nrow(M)+1))
for(i in 1:ncol(M)){
  for(j in 1:nrow(M)){
    points(i,j,col=1,pch=circle[M[j,i]+1],cex=1)
  }
}

# plot horizontal line
for(j in 1:nrow(M)){    
  abline(h=j)    
}



# scenarioF MHL=0.11, AML= 0.5
hapinfo<-c(c(rep("TTTTTT",90),"CTTTT","CTTTC","TCTTC","TTTTC","TCTTT","TTTCT","TTTC","TCT","T"),rep("CCCCCC",1))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
write.table(data.frame(hapinfo),file="6.txt",quote=F)


