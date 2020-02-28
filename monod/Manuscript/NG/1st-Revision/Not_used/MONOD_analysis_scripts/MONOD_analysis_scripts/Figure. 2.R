#!/usr/bin/R
# Methylation Haploinfo to Methylation haplotype load (mode weigth i)
# Contact: Shicheng Guo
# Version 1.3
# Update: May/29/2015

hap_comb_tmp<-function(hapinfo.single){
  # hap_com.single is a character of haplotype, such as "TCTCTT"
  hap_com.single<-c()
  for(i in 1:nchar(hapinfo.single)){
    for(j in i:(nchar(hapinfo.single))){
      hap_com.single<-c(hap_com.single,(substr(hapinfo.single,i,j)))
    } 
  }
  return(hap_com.single)
}

hap_comb<-function(hapinfo){
  # haplotype is array of observed methylation haplotype and return
  hap_comb<-unlist(lapply(hapinfo,hap_comb_tmp))
  return(hap_comb)
}

mhl<-function(hap_comb_array){
  mhl_value<-c()
  hap_comb_array<-unlist(hap_comb_array)
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  unmeth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="T"
  hap_comb_array<-c(hap_comb_array[meth_hap],hap_comb_array[unmeth_hap])
  for(i in unique(nchar(hap_comb_array))){
    input_tmp<-which(nchar(hap_comb_array)==i)  
    nmeth<-length(grep("C",hap_comb_array[input_tmp]))
    nunmeth<-length(grep("T",hap_comb_array[input_tmp]))
    mhl_value<-c(mhl_value,i*nmeth/(nunmeth+nmeth))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}

mf<-function(hapinfo){
  mf_value<-c()
  hap_comb_array<-unlist(hapinfo)
  meth_hap=unlist(lapply(hap_comb_array,function(x) unlist(strsplit(x,""))))
  nmeth<-length(grep("C",meth_hap))
  nunmeth<-length(grep("T",meth_hap))
  mf_value<-nmeth/(nunmeth+nmeth)
  return(mf_value)
}

# Figuure 2A
hapinfo<-c(rep("TTTT",16))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2B
hapinfo<-c(rep("CCCC",16))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2C
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
MHL<-mhl(hap_comb(unique(hapinfo)))
MF<-mf(hapinfo)
me<-methentropy(hapinfo)
# Figuure 2D
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
MHL<-mhl(hap_comb(unique(hapinfo)))
# Figuure 2E
hapinfo<-c(rep("CCTT",8),rep("TTCC",8))
MHL<-mhl(hap_comb(unique(hapinfo)))

