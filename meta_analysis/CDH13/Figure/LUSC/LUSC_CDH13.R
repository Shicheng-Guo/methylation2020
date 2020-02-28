setwd("/media//puweilin/新加卷/Gengxin/")
anno=read.csv("GPL13534_HumanMethylation450_15017482_v.1.1.csv",header=T,sep="\t")
#data preparation
setwd("/media//puweilin/新加卷/Gengxin/LUSC/")
load("LUSC_CDH13.RData")
LUSC_CDH13=LUSC_CDH13[,which(LUSC_CDH13[4,] %in% c("01","11"))]
type=LUSC_CDH13[4,]
for(i in 1:length(type)) {if(type[i]=="01") {type[i]="Case"} else if(type[i]=="11") {type[i]="Control"}  }
LUSC_CDH13=na.omit(LUSC_CDH13[-c(1:4),])
LUSC_CDH13=t(LUSC_CDH13)
LUSC_CDH13=apply(LUSC_CDH13,2,function(x) as.numeric(x))
seq=match(colnames(LUSC_CDH13),anno[,1])
cg_ref=anno[seq,c(1,15,23,24,25)]
cg_ref2=cg_ref[order(-as.numeric(cg_ref[,2])),]  # cg_ref2 is the reference information for the cgsites of CDH13
cg_ref_CDH13=cg_ref2  
save(cg_ref_CDH13,file="cg_ref_CDH13.RData")
seq=match(cg_ref2[,1],colnames(LUSC_CDH13))
new_LUSC_CDH13=LUSC_CDH13[,seq]
type_total=rep(type,dim(new_LUSC_CDH13)[2])
methy_total=c();for(i in 1:dim(new_LUSC_CDH13)[2]){ methy_total=append(methy_total,new_LUSC_CDH13[,i])}
cgsite_total=c(); for(i in 1:dim(new_LUSC_CDH13)[2]) { cgsite_total=append(cgsite_total,rep(colnames(new_LUSC_CDH13)[i],dim(new_LUSC_CDH13)[1]))}
LUSC_CDH13_newdata=as.data.frame(cbind(type_total,cgsite_total,methy_total))
LUSC_CDH13_newdata[,3]=as.numeric(as.character(LUSC_CDH13_newdata[,3]))  # LUSC_CDH13_newdata is the data processed to use for the data summarization 
rownames(LUSC_CDH13_newdata)=c()
mean=aggregate(LUSC_CDH13_newdata$methy_total,by=list(type=LUSC_CDH13_newdata$type_total,cgsite=LUSC_CDH13_newdata$cgsite_total),mean)
CI1=aggregate(LUSC_CDH13_newdata$methy_total,by=list(type=LUSC_CDH13_newdata$type_total,cgsite=LUSC_CDH13_newdata$cgsite_total),function(x) {t.test(x)$conf[1]})
CI2=aggregate(LUSC_CDH13_newdata$methy_total,by=list(type=LUSC_CDH13_newdata$type_total,cgsite=LUSC_CDH13_newdata$cgsite_total),function(x) {t.test(x)$conf[2]})
LUSC_CDH13_forplot=data.frame(cbind(mean,CI1[,3],CI2[,3]))     # LUSC_CDH13_forplot is the dataframe which summarizes the mean and confidence interval of CDH13 methylation cgsites 
colnames(LUSC_CDH13_forplot)=c("type","cgsite","mean","CI1","CI2")  
save(LUSC_CDH13_forplot,file="LUSC_CDH13_forplot.RData")

#LUSC_total data for ploting
LUSC_CDH13_forplot$cgsite=factor(LUSC_CDH13_forplot$cgsite,order=T,levels=cg_ref2[,1]) #order the cgsites to be the same as in the chromesome  
#ggplot

ggplot(LUSC_CDH13_forplot, aes(x=cgsite, y=mean, colour=type, group=type,shape=type)) + 
  geom_errorbar(aes(ymin=CI1, ymax=CI2), colour="black", width=.4,lwd=0.8) +
  geom_line(lwd=0.8)+
  geom_point( size=4, fill="white")+ # 21 is filled circle
  xlab("CpG site")+
  ylab("DNA methylation(β value)")+scale_y_continuous(limits=c(0, 1),    # Set y range
                                                      breaks=c(0,0.2,0.4,0.6,0.8,1)) +                       # Set tick every 0.2
  theme_bw() +
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





#LUSC_island data for ploting 
island_cg=cg_ref2[which(cg_ref2[,5] %in% c("N_Shore","Island","S_Shore")),1]
LUSC_CDH13_forplot_island=LUSC_CDH13_forplot[which(LUSC_CDH13_forplot$cgsite %in% island_cg),]
LUSC_CDH13_forplot_island$cgsite=factor(LUSC_CDH13_forplot_island$cgsite,order=T,levels=island_cg)


#ggplot ploting 
ggplot(LUSC_CDH13_forplot_island, aes(x=cgsite, y=mean, colour=type, group=type,shape=type)) + 
  geom_errorbar(aes(ymin=CI1, ymax=CI2), colour="black", width=.4,lwd=0.8) +
  geom_line(lwd=0.8) +
  geom_point( size=4, fill="white") + # 21 is filled circle
  xlab("CpG site") +
  ylab("DNA methylation(β value)")  +
  scale_y_continuous(limits=c(0, 1),    # Set y range
                     breaks=c(0,0.2,0.4,0.6,0.8,1)) +                       # Set tick every 0.2
  theme_bw() +
  theme( legend.position=c(0.1,0.9))+
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





# LUSC_CDH13_nonisland_unmatched
non_island_cg=cg_ref2[-which(cg_ref2[,5] %in% c("N_Shore","Island","S_Shore")),1]
LUSC_CDH13_forplot_non_island=LUSC_CDH13_forplot[which(LUSC_CDH13_forplot$cgsite %in% non_island_cg),]
LUSC_CDH13_forplot_non_island$cgsite=factor(LUSC_CDH13_forplot_non_island$cgsite,order=T,levels=non_island_cg)

#ggplot ploting 
ggplot(LUSC_CDH13_forplot_non_island, aes(x=cgsite, y=mean, colour=type, group=type,shape=type)) + 
  geom_errorbar(aes(ymin=CI1, ymax=CI2), colour="black", width=.4,lwd=0.8) +
  geom_line(lwd=0.8)+
  geom_point( size=4, fill="white")+ # 21 is filled circle
  xlab("CpG site")+
  ylab("DNA methylation(β value)")+scale_y_continuous(limits=c(0, 1),    # Set y range
                                                      breaks=c(0,0.2,0.4,0.6,0.8,1)) + # Set tick every 0.2
  theme_bw() +
  theme( legend.position=c(0.1,0.1))+
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


