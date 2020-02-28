setwd("/home/wangjc01/wlpu/methy_RNA/ESCA/ESCA/ESCA_methylation")
ESCA_methy=read.table("ESCA_methylation.txt",header=T,sep="\t")
ESCA_methy=ESCA_methy[-1,]
cgsite=ESCA_methy[,1]
rownames(ESCA_methy)=cgsite
ESCA_methy=ESCA_methy[,-1]
temp=as.character(colnames(ESCA_methy))
barcode1=as.character(sapply(temp,function(x) strsplit(x,"[.]")[[1]][1]))
barcode2=as.character(sapply(temp,function(x) strsplit(x,"[.]")[[1]][2]))
barcode3=as.character(sapply(temp,function(x) strsplit(x,"[.]")[[1]][3]))
barcode=paste(barcode1,barcode2,barcode3,sep="-")
type1=as.character(sapply(temp,function(x) strsplit(x,"[.]")[[1]][4]))
type=as.character(sapply(type1,function(x) substr(x,1,2)))
ESCA_methy=rbind(type,as.matrix(ESCA_methy))
ESCA_methy=na.omit(ESCA_methy)
colnames(ESCA_methy)=barcode
seq=ESCA_methy[1,] %in% c("01","11")
ESCA_methy=ESCA_methy[,seq]
save(ESCA_methy,file="ESCA_methy.RData")
barcode_times=as.data.frame(table(colnames(ESCA_methy)))
barcode_times=subset(barcode_times,barcode_times[,2]==2);
unique_barcode =unique(as.vector(barcode_times[,1]));
seq=which(colnames(ESCA_methy) %in% unique_barcode);
ESCA_methy_matched=ESCA_methy[,seq]
ESCA_methy_matched = ESCA_methy_matched[-381166,];

methy_ttest=function(data,FDR=TRUE,cutoff=0.3){
  type01=which(data[1,]=="01")
  type11=which(data[1,]=="11")
  test.result=apply(data[-1,],1,function(x) {a=t.test(as.numeric(x[type01]),as.numeric(x[type11])); b=as.numeric(a$p.value); c=as.numeric(data.frame(a$estimate)[1:2,1]);  e=ifelse(as.numeric(c[1])>as.numeric(c[2]),"Hyper","Hypo"); d=append(append(b,c),e); return(d)  })
  rownames(test.result)=c("pvalue","cancer_avg","normal_avg","HyperorHypo")
  test.result=t(test.result)
  test.result=as.data.frame(test.result)
  rownames(test.result)=rownames(data)[-1]
  test.result[,1:3]=apply(test.result[,1:3],2,function(x) as.numeric(x))
  if(FDR){
    p.adjusted = p.adjust(test.result[,1],method = "fdr");
    test.result = data.frame(test.result,p.adjusted);
    test.selected=subset(test.result,test.result[,5]<0.05 & abs(test.result[,2]-test.result[,3])>cutoff)
    test.selected=test.selected[order(test.selected[,1]),]
    return(test.selected)
    
  }else{
    test.selected=subset(test.result,test.result[,1]<(0.05/dim(test.result)[1]) & abs(test.result[,2]-test.result[,3])>cutoff)
    test.selected=test.selected[order(test.selected[,1]),]
    return(test.selected)
  }
}


ESCA_methytest=methy_ttest(ESCA_methy_matched)
save(ESCA_methytest,file="ESCA_methytest.RData")


ESCA_testdata=ESCA_methy_matched[c(1,match(rownames(ESCA_methytest),rownames(ESCA_methy_matched))),]
ESCA_testdata=as.data.frame(t(ESCA_testdata))
ESCA_testdata[,1]=as.factor(ESCA_testdata[,1])
ESCA_testdata[,-1]=apply(ESCA_testdata[,-1],2,function(x) as.numeric(x))
ESCA_testdata=subset(ESCA_testdata,as.character(ESCA_testdata[,1]) %in% c("01","11"))
ESCA_testdata[,1]=as.factor(as.character(ESCA_testdata[,1]))
save(ESCA_testdata,file="ESCA_testdata.RData")

library(randomForest)
rf_ESCA=randomForest(ESCA_testdata[,-1],ESCA_testdata[,1],proximity=T,importance=T,ntree=1500)
save(rf_ESCA,file="rf_ESCA.RData")
ESCA_imp=rf_ESCA$importance
ESCA_imp=ESCA_imp[order(-ESCA_imp[,3]),]
cgsite_top10=rownames(ESCA_imp)[1:10]
ESCA_top10=ESCA_testdata[,c(1,match(cgsite_top10,colnames(ESCA_testdata)))]
rf_ESCA_top10=randomForest(ESCA_top10[,-1],ESCA_top10[,1],proximity=T,importance=T,ntree=1500)
save(rf_ESCA_top10,file="rf_ESCA_top10.RData");
annotation=read.csv("../GPL13534_HumanMethylation450_15017482_v.1.1.csv",header=T,sep="\t")


cgsite=rownames(ESCA_imp)
seq=match(cgsite,annotation[,1])
cgsite=annotation[seq,c(1,23,25)]
examplefile=read.table("../Raw data/merthylation/DNA_Methylation/JHU_USC__HumanMethylation450//Level_3//jhu-usc.edu_ESCA.HumanMethylation450.10.lvl-3.TCGA-V5-A7RC-01B-11D-A409-05.txt",header=T,sep="\t")
gene=examplefile[match(cgsite[,1],examplefile[,1]),3]
cgsite=cbind(cgsite,gene)
save(cgsite,file="cgsite.RData")


library(xlsx)
cg_top2000 = cgsite[1:2000,];
genes = as.data.frame(table(cgsite_top2000[,4]))
genes = genes[order(-genes[,2]),]
genes = subset(genes, genes[,2]>=4);
write.xlsx2(genes,file="Significant_Genes.xlsx",sheetName = "Significant Genes")




#for plotting the gene distribution between case and control 
library(ggplot2)
get_gene_data=function(x){
  cg=examplefile[which(examplefile[,3]==x),1]
  seq=match(cg,rownames(ESCA_methy_matched))
  data=(ESCA_methy_matched[seq,])
  data=t(na.omit(data))
  cg=colnames(data)
  pos=annotation[match(cg,annotation[,1]),c(15,24,25)]
  cg=cbind(cg,pos)   
  cgsite_ordered=cg[order(cg[,2]),1] #return1
  pos=cg[order(cg[,2]),3] #return2
  cg=cg[order(cg[,2]),] #return3
  type=ESCA_methy_matched[1,]
  type_total=rep(type,dim(data)[2])
  cgsite_total=c(); for(i in 1:dim(data)[2]) { cgsite_total=append(cgsite_total,rep(colnames(data)[i],length(type)))}
  methy_total=c(); for (i in 1:dim(data)[2]) { methy_total=append(methy_total,data[,i])}
  data_new=data.frame(type_total,cgsite_total,methy_total)
  data_new[,3]=as.numeric(as.character(data_new[,3]))
  mean=aggregate(data_new$methy_total,by=list(type=data_new$type_total,cgsite=data_new$cgsite_total),mean)
  CI1=aggregate(data_new$methy_total,by=list(type=data_new$type_total,cgsite=data_new$cgsite_total),function(x) {t.test(x)$conf[1]})
  CI2=aggregate(data_new$methy_total,by=list(type=data_new$type_total,cgsite=data_new$cgsite_total),function(x) {t.test(x)$conf[2]})
  data_forplot=data.frame(cbind(mean,CI1[,3],CI2[,3]))
  colnames(data_forplot)=c("type","cgsite","mean","CI1","CI2") 
  data_forplot$cgsite=factor(data_forplot$cgsite,order=T,levels=cgsite_ordered)#return4
  library(ggplot2)
  p=theme_bw() +
    theme( legend.position=c(0.05,0.1))+
    theme(legend.title = element_text(colour="black", size=20, face="bold"))+
    theme(legend.text = element_text(colour="black", size = 16, face = "bold"))+
    theme(panel.grid.major=element_line(colour=NA))+ 
    theme(plot.background = element_blank(),panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank(),panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black'))+
    theme(axis.text.x=element_text(size=15,vjust=0.9,face="bold",hjust=1,angle=45))+
    theme(axis.title.x=element_text(size=25,vjust=-0.1))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.y=element_text(size=22,vjust=0.3))#draws x and y axis line
  q=ggplot(data_forplot, aes(x=cgsite, y=mean, colour=type, group=type,shape=type)) + 
    geom_errorbar(aes(ymin=CI1, ymax=CI2), colour="black", width=.4,lwd=0.8) +geom_line(lwd=0.8)+
    geom_point( size=4, fill="white")+ # 21 is filled circle
    xlab("CpG site")+
    ylab("DNA methylation(beta value)")+scale_y_continuous(limits=c(0, 1),                  # Set y range
                                                           breaks=c(0,0.2,0.4,0.6,0.8,1)) +p                 # Set y range
  final=list(data,cg,data_forplot,q)
  return(final)
}


# plotting the significant genes previously saved
sig_genes = read.xlsx2("Significant_Genes.xlsx",sheetIndex = 1,stringsAsFactors=F);
for( gene in sig_genes[,1]){
  Gene_Data = get_gene_data(gene);
  save(Gene_Data,file=paste(gene,".Rdata",sep=""));
  ggsave(paste(gene,".jpeg",sep=""),Gene_Data[[4]],width = 20, height = 10)
}


#get the RNA-seq dataset from TCGA and find out if some of the previous significant genes' methylation status are associated with the expression
load("ESCA_rnaseq.RData");
library(ggplot2)
q=theme_bw() +
  theme(panel.grid.major=element_line(colour=NA))+ 
  theme(plot.background = element_blank(),panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank(),panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))+
  theme(axis.text.x=element_text(size=15,vjust=0.9,face="bold"))+
  theme(axis.title.x=element_text(size=25,vjust=-0.1))+
  theme(axis.text.y=element_text(size=20))+
  theme(axis.title.y=element_text(size=22,vjust=0.3))+
  theme(legend.position=c(0.95,0.5))+
  theme(legend.title = element_text(size=20, face="bold"))+
  theme(legend.text = element_text(size = 16, face = "bold"))


expression_diff=function(gene_name,rnaseq_data){
  gene_list = rownames(rnaseq_data)[-1];
  seq=match(gene_name,gene_list)
  if(!is.na(seq)){
  gene_data=rnaseq_data[c(1,seq),]
  seq=which(gene_data[1,] %in% c("01","11"))
  gene_data=gene_data[,seq]
  gene_data=t(gene_data)
  gene_rownames=as.data.frame(table(rownames(gene_data)))
  gene_rownames=gene_rownames[order(-gene_rownames[,2]),]
  gene_rownames=subset(gene_rownames,as.numeric(gene_rownames[,2])==2)
  seq=which(rownames(gene_data) %in% as.character(gene_rownames[,1]))
  gene_data=gene_data[seq,]
  gene_rownames=rownames(gene_data)
  gene_rownames1=paste(rownames(gene_data),gene_data[,1],sep="-")
  gene_data=data.frame(type=gene_data[,1],gene_expression=as.numeric(gene_data[,2]))
  rownames(gene_data)=gene_rownames1
  case_seq=which(as.character(gene_data[,1])=="01")
  control_seq=which(as.character(gene_data[,1])=="11")
  data_case=gene_data[case_seq,2]
  data_control=gene_data[control_seq,2]
  pvalue=wilcox.test(data_case,data_control)$p.value
  plot=ggplot(gene_data)+geom_jitter(aes(type,gene_expression,colour=type))+q
  result=list(gene_data,pvalue,plot)
  savename1=paste(gene_name,"data.csv",sep="_")
  savename2=paste(gene_name,"pvalue.txt",sep="_")
  savename3=paste(gene_name,".jpeg",sep="")
  write.csv(gene_data,file=savename1,quote=F)
  write.table(pvalue,file=savename2,quote=F)
  ggsave(file=savename3,plot=plot)
  return(result)}
}


for( gene in sig_genes[,1]){
  Expression_data = expression_diff(gene,ESCA_rnaseq);
  save(Expression_data, file=paste(gene,".Rdata",sep=""))
}



# we got the GEO dataset which consists of 4 cases and 8 controls and we want to do validation for the TCGA methylation dataset 

setwd("f:/GEO datasets/ESCA/");
GEO =read.csv("GSE52826_series_matrix.csv",header=T)
rownames(GEO) =GEO[,1];
GEO =GEO[,-1];
clinical = read.csv("clinical.csv")
setwd("f:/Esphageal Cancer/2016.2.26/")

library(ggplot2)
get_gene_data=function(x){
  cg=examplefile[which(examplefile[,3]==x),1]
  seq=match(cg,rownames(GEO));
  data=(GEO[seq,])
  data=t(na.omit(data))
  cg=colnames(data)
  pos=annotation[match(cg,annotation[,1]),c(15,24,25)]
  cg=cbind(cg,pos)   
  cgsite_ordered=cg[order(cg[,2]),1] #return1
  pos=cg[order(cg[,2]),3] #return2
  cg=cg[order(cg[,2]),] #return3
  type=clinical[,2]
  type_total=rep(type,dim(data)[2])
  cgsite_total=c(); for(i in 1:dim(data)[2]) { cgsite_total=append(cgsite_total,rep(colnames(data)[i],length(type)))}
  methy_total=c(); for (i in 1:dim(data)[2]) { methy_total=append(methy_total,data[,i])}
  data_new=data.frame(type_total,cgsite_total,methy_total)
  data_new[,3]=as.numeric(as.character(data_new[,3]))
  mean=aggregate(data_new$methy_total,by=list(type=data_new$type_total,cgsite=data_new$cgsite_total),mean)
  CI1=aggregate(data_new$methy_total,by=list(type=data_new$type_total,cgsite=data_new$cgsite_total),function(x) {t.test(x)$conf[1]})
  CI2=aggregate(data_new$methy_total,by=list(type=data_new$type_total,cgsite=data_new$cgsite_total),function(x) {t.test(x)$conf[2]})
  data_forplot=data.frame(cbind(mean,CI1[,3],CI2[,3]))
  colnames(data_forplot)=c("type","cgsite","mean","CI1","CI2") 
  data_forplot$cgsite=factor(data_forplot$cgsite,order=T,levels=cgsite_ordered)#return4
  library(ggplot2)
  p=theme_bw() +
    theme( legend.position=c(0.05,0.1))+
    theme(legend.title = element_text(colour="black", size=20, face="bold"))+
    theme(legend.text = element_text(colour="black", size = 16, face = "bold"))+
    theme(panel.grid.major=element_line(colour=NA))+ 
    theme(plot.background = element_blank(),panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank(),panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black'))+
    theme(axis.text.x=element_text(size=15,vjust=0.9,face="bold",hjust=1,angle=45))+
    theme(axis.title.x=element_text(size=25,vjust=-0.1))+
    theme(axis.text.y=element_text(size=20))+
    theme(axis.title.y=element_text(size=22,vjust=0.3))#draws x and y axis line
  q=ggplot(data_forplot, aes(x=cgsite, y=mean, colour=type, group=type,shape=type)) + 
    geom_errorbar(aes(ymin=CI1, ymax=CI2), colour="black", width=.4,lwd=0.8) +geom_line(lwd=0.8)+
    geom_point( size=4, fill="white")+ # 21 is filled circle
    xlab("CpG site")+
    ylab("DNA methylation(beta value)")+scale_y_continuous(limits=c(0, 1),                  # Set y range
                                                           breaks=c(0,0.2,0.4,0.6,0.8,1)) +p                 # Set y range
  final=list(data,cg,data_forplot,q)
  return(final)
}


# plotting the significant genes previously saved
sig_genes = read.xlsx2("Significant_Genes.xlsx",sheetIndex = 1,stringsAsFactors=F);
for( gene in sig_genes[,1]){
  Gene_Data = get_gene_data(gene);
  save(Gene_Data,file=paste(gene,".Rdata",sep=""));
  ggsave(paste(gene,".jpeg",sep=""),Gene_Data[[4]],width = 20, height = 10)
}

