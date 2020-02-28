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

###############################################################################################################################
### methylation haplotype load for the Different situation/haplotypes (6 situation) AML: absolute methylation level (average)
###############################################################################################################################
# scenarioA MHL=1, AML= 1     # mhl:(1/10+8/60)
hapinfo<-c(rep("CCCC",8),rep("CCCC",8))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
write.table(data.frame(hapinfo),file="1.txt",quote=F,col.names=F,row.names=F)
# scenarioB MHL=0, AML= 0     # mhl:(1/10+8/60)
hapinfo<-c(rep("TTTT",8),rep("TTTT",8))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
write.table(data.frame(hapinfo),file="2.txt",quote=F,col.names=F,row.names=F)

# scenarioC MHL=0.5, AML= 0.5
hapinfo<-c(rep("CCCC",8),rep("TTTT",8))
mhl_rlt<-mhl(hap_comb((hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
write.table(data.frame(hapinfo),file="3.txt",quote=F,col.names=F,row.names=F)

# scenarioD MHL=0.1625, AML= 0.5
hapinfo<-c("TCCT","TTTC","CTTC","CTCT" ,"TTCC", "CTTT", "CCTT" ,"TTCT" ,"CCCC" ,"CCTC" ,"TCTT" ,"TTTT" ,"CTCC", "TCTC", "TCCC" ,"CCCT")
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
write.table(data.frame(hapinfo),file="4.txt",quote=F,col.names=F,row.names=F)

# scenarioE MHL=0.23, AML= 0.5     # mhl:(1/10+8/60)
hapinfo<-c(rep("CCTT",8),rep("TTCC",8))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
write.table(data.frame(hapinfo),file="5.txt",quote=F,col.names=F,row.names=F)

# scenarioF MHL=0.11, AML= 0.5
hapinfo<-c(rep("CCTT",8),rep("CCTT",8))
mhl_rlt<-mhl(hap_comb(unique(hapinfo)))
aml_rlt<-hap2methfreq(hapinfo)
mhl_rlt
aml_rlt
write.table(data.frame(hapinfo),file="6.txt",quote=F,col.names=F,row.names=F)


# methylation circle plot 
hapinfo_matrix <- data.frame(CpG=rep(1:4,16), read=rep(1:16, each=4), value=c(rep(1,32),rep(16,32)))
library("ggplot2")
ggplot() + 
  scale_shape_identity() + 
  geom_point(data=hapinfo_matrix, aes(x=CpG, y=read, shape=value), size=7)+
  coord_cartesian(ylim = c(0, 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


