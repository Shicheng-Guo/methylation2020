#######################################################################################################################
###   Title : Methylation haplotype load and absolute methylation level
###   Lab:  Kun's Lab (http://zhang.openwetware.org/)   
###   Time :  Sep/29/2015 
###   Update: re-define MHL to identify continous long methylation and non-methylation reads
#######################################################################################################################

hap2matrix<-function(hapinfo){
  hapinfo.array<-lapply((hapinfo),function(x) unlist(strsplit(x,"")))
  y<-unlist(hapinfo.array)
  hapinfo.array<-matrix(y,ncol=length(x[[1]]),byrow=T)
  return(hapinfo.array)
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
    mhl_value<-c(mhl_value,i*(nmeth+nunmeth)/sum(nchar(hap_comb_array)==i))
    print(sum(nchar(hap_comb_array)==i))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}


# MHL=0.325, AML= 0.5
hapinfo<-c("TCCT","TTTC","CTTC","CTCT" ,"TTCC", "CTTT", "CCTT" ,"TTCT" ,"CCCC" ,"CCTC" ,"TCTT" ,"TTTT" ,"CTCC", "TCTC", "TCCC" ,"CCCT")
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt

matrix<-data.matrix(data.frame(hap2matrix(hapinfo)))
matrix<-matrix-1
cor<-cor(t(matrix))
cor[lower.tri(cor,diag=T)]<-NA
cor<-na.omit(as.vector(cor))
mean(cor)

# MHL=0.23 (8/60), AML= 0.5
hapinfo<-c(rep("CCTT",8),rep("TTCC",8))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt

# MHL=0.23, AML= 0.5
hapinfo<-c(rep("CCTT",8),rep("CCTT",8))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt

# MHL=0.5, AML= 0.5
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
mhl_rlt<-mhl(hap_comb((hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
# MHL=1, AML= 0.5
hapinfo<-c(rep("CCCC",16),rep("TTTT",0))
mhl_rlt<-mhl(hap_comb((hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
# MHL=1, AML= 0.5
hapinfo<-c(rep("CCCC",0),rep("TTTT",16))
mhl_rlt<-mhl(hap_comb((hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt








hapinfo_matrix<-as.matrix(unique(hapinfo),ncol=1)
hapinfo_matrix<-as.matrix(hapinfo,ncol=1)
hapinfo_matrix




mhl_rlt<-mhl(hap_comb(hapinfo))
mhl_rlt



haplo2ggplot<-function(x){
  
  
}


hapinfo_matrix <- data.frame(CpG=rep(1:4,16), read=rep(1:16, each=4), value=c(rep(1,32),rep(16,32)))
library("ggplot2")
ggplot() + 
  scale_shape_identity() + 
  geom_point(data=hapinfo_matrix, aes(x=CpG, y=read, shape=value), size=7)+
  coord_cartesian(ylim = c(0, 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))





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
  # haplotype is array of observed methylation haplotype
  hap_comb<-unlist(lapply(hapinfo,hap_comb_tmp))
  return(hap_comb)
}

mhl<-function(hap_comb_array){
  # hap_comb_array is combinarial methylation haplotype
  mhl_value<-c()
  hap_comb_array<-unlist(hap_comb_array)
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  unmeth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="T"
  hap_comb_array<-c(hap_comb_array[meth_hap],hap_comb_array[unmeth_hap])
  for(i in unique(nchar(hap_comb_array))){
    input_tmp<-which(nchar(hap_comb_array)==i)  
    nmeth<-length(grep("C",hap_comb_array[input_tmp]))
    nunmeth<-length(grep("T",hap_comb_array[input_tmp]))
    mhl_value<-c(mhl_value,i*(nmeth+numeth)/length(hap_comb_array))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}



hap2methfreq<-function(hapinfo){
  hapinfo<-paste(hapinfo,collapse="")
  hapinfo<-strsplit(hapinfo,"")
  num_C<-length(grep("C",unlist(hapinfo)))
  num_T<-length(grep("T",unlist(hapinfo)))
  methfreq<-  num_C/(num_C+num_T)
  return(methfreq)
}

# MHL=0.5
hapinfo<-c()
for(i in 1:500){
  hapinfo<-c(hapinfo,paste(sample(c("C","T"),4,replace=T),collapse = ""))
}

mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
mhl_rlt
aml_rlt<-hap2methfreq(hapinfo)
aml_rlt

# eventally, I find a classic example to show mhl

hapinfo<-c("TTCC","CCCC","CTCT")


# MHL=0.5
hapinfo<-c(rep("CCTT",8),rep("TTCC",9))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
# MHL=0.5
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
mhl_rlt

# MHL=1
hapinfo<-c(rep("CCCC",16),rep("TTTT",0))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
mhl_rlt

# MHL=0
hapinfo<-c(rep("CCCC",0),rep("TTTT",16))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
mhl_rlt


hapinfo_matrix<-as.matrix(unique(hapinfo),ncol=1)


hapinfo_matrix<-as.matrix(hapinfo,ncol=1)
hapinfo_matrix




mhl_rlt<-mhl(hap_comb(hapinfo))
mhl_rlt



haplo2ggplot<-function(x){
  
  
}


hapinfo_matrix <- data.frame(CpG=rep(1:4,16), read=rep(1:16, each=4), value=c(rep(1,32),rep(16,32)))
library("ggplot2")
ggplot() + 
  scale_shape_identity() + 
  geom_point(data=hapinfo_matrix, aes(x=CpG, y=read, shape=value), size=7)+
  coord_cartesian(ylim = c(0, 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))





source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
library("Rsamtools")

  
  
  
  
  


