setwd("f:/Esphageal Cancer/2016.3.17/")
load("ESCA_methy_matched.RData")
source("site2region.R")
source("testPvalue.R")
ESCA_methytest = methy_ttest(ESCA_methy_matched)
save(ESCA_methytest,file="ESCA_methytest.RData")

ESCA_group = sites2region(ESCA_methytest,window = 6)
save(ESCA_group, file = "ESCA_group.RData")

Region_stats = sapply(ESCA_group, FUN = function(x){
  Pvalue = combineTestPvalsMeth(ESCA_methytest[x,1])
  McaM = mean(ESCA_methytest[x,2])
  McoM = mean(ESCA_methytest[x,3])
  MeanDiff = McaM - McoM
  return(c(McaM,McoM, MeanDiff, Pvalue))
})

Region_stats = as.data.frame(t(Region_stats))
FDR = p.adjust(Region_stats[,4],method = "bonferroni")
Region_stats = data.frame(Region_stats,FDR)
colnames(Region_stats) = c("McaM", "McoM", "MeanDiff", "Pvalue","FDR")
save(Region_stats, file="Region_stats.RData")

#select the cutoff for the region selection
seq = which( Region_stats$FDR <=0.01 & Region_stats$MeanDiff >0.40 & Region_stats$McoM <0.30)
Regions_selected = Region_stats[seq,]
Regions_selected = Regions_selected[order(Regions_selected$FDR),]

# match the regions to the genes 
examplefile=read.table("../Raw data/merthylation/DNA_Methylation/JHU_USC__HumanMethylation450//Level_3/jhu-usc.edu_ESCA.HumanMethylation450.10.lvl-3.TCGA-V5-A7RC-01B-11D-A409-05.txt",skip =1, header=T,sep="\t")
genes = sapply( rownames(Regions_selected), FUN =function(x){
  seq = which(names(ESCA_group) == x)
  cgsite_seq = ESCA_group[[seq]]
  cgsite = rownames(ESCA_methytest)[cgsite_seq]
  genes = examplefile[match(cgsite,examplefile[,1]),3]
  gene = unique(genes)
  gene_combine = paste(gene, collapse="&")
  return(gene_combine)
})
Gene = as.vector(unlist(genes))
Regions_selected = data.frame(Regions_selected, Gene)

Location = sapply( rownames(Regions_selected), FUN =function(x){
  seq = which(names(ESCA_group) == x)
  cgsite_seq = ESCA_group[[seq]]
  cgsite = rownames(ESCA_methytest)[cgsite_seq]
  cgsite_first = cgsite[1]
  cgsite_last = cgsite[length(cgsite)]
  locus_first = examplefile[match(cgsite_first,examplefile[,1]),5]
  locus_last = examplefile[match(cgsite_last,examplefile[,1]),5]
  Locus = paste(locus_first,locus_last,sep=" - ")
  return(Locus)
})
Locus = as.vector(unlist(Location))
CHR = sapply(rownames(Regions_selected), function(x) { return(strsplit(as.character(x),"_")[[1]][1])})
CHR = as.vector(CHR)
Regions_selected = data.frame(Regions_selected,CHR,Locus)
save(Regions_selected, file="Regions_selected.RData")



#Regional plotting 
library(ggplot2)
get_region_data=function(region){
  seq = which(names(ESCA_group) == region)
  cgsite_seq = ESCA_group[[seq]]
  cgsite = rownames(ESCA_methytest)[cgsite_seq]
  seq=match(cgsite,rownames(ESCA_methy_matched))
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
    theme( legend.position=c(0.9,0.9))+
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
for( region in rownames(Regions_selected)){
  Region_Data = get_region_data(region);
  ggsave(paste(region,".jpeg",sep=""),Region_Data[[4]],width = 20, height = 10)
  seq = which(names(ESCA_group) == region)
  cgsite_seq = ESCA_group[[seq]]
  data = ESCA_methytest[cgsite_seq,]
  write.xlsx2(data,file= paste(region,".xlsx",sep=""),sheetName = "region")
}

