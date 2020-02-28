########################################################################################
###   Title: Group Specificity Index (GSI) for Genome-wide Methylation Haplotype dataset
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 9/1/2015
########################################################################################

library("impute")

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
} 
bedwithgap<-function(bed,gap){
  bed<-as.matrix(bed)
  bed[,2]=as.numeric(bed[,2])-gap
  bed[,3]=as.numeric(bed[,3])+gap
  bed<-data.frame(bed)
  bed
}

Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(as.character(cor),split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}

bed2cor<-function(bed){
  cor<-apply(bed,1,function(x){paste(unlist(strsplit(x,"\t"))[1],":",unlist(strsplit(x,"\t"))[2],"-",unlist(strsplit(x,"\t"))[3],sep="")})
  cor<-gsub("[ ]","",cor)
  return(cor)
}


#########################################################################################

setwd("/home/shg047/monod/oct/data")
file1<-read.table("WGBS_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
file1<-data.matrix(file1[,13:56])  # remove H1 and WBC, cancer
colnames(file1)
colnames(file1)<-gsub("_","-",colnames(file1))
colname2<-unlist(lapply(colnames(file1),function(x) unlist(strsplit(x,"[.]"))[1]))
saminfo2<-read.table("/home/shg047/monod/phase2/newsaminfo.txt",head=T,sep="\t",as.is=T)
saminfo2<-saminfo2[na.omit(match(colname2,saminfo2[,1])),]
colnames(file1)<-saminfo2[,2]
colnames(file1)

f2<-RawNARemove(file1,missratio=0.4)
file1<-impute.knn(f2)$data

group=names(table(colnames(file1)))
index=colnames(file1)
gsi<-c()
gmaxgroup<-c()
for(i in 1:nrow(file1)){
gsit<-0
gmax<-names(which.max(tapply(as.numeric(file1[i,]),index,mean)))
for(j in 1:length(group)){
  tmp<-(1-10^(mean(file1[i,][which(index==group[j])]))/10^(mean(file1[i,][which(index==gmax)])))/(length(group)-1)
  gsit<-gsit+tmp
}
gmaxgroup<-c(gmaxgroup,gmax)
gsi<-c(gsi,gsit)
print(c(gmax,gsit))
}
rlt=data.frame(region=rownames(file1),group=gmaxgroup,GSI=gsi)
write.table(rlt,file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",col.names=T,row.names=F,quote=F,sep="\t")

# each take top 5 tissue-specific methylation regions.
data<-read.table(file="Table.GSI.WGBS.Remove.H1.WBC.rlt.txt",head=T,sep="\t",as.is=T)
# heat only have 3 high GSI regions
sum(table(subset(data,GSI>0.5)[,2]))
head(data)
tissue<-names(table(data[,2]))
tissue<-sort(c("Brain","Heart","muscle","Vessel","Spleen","Kidney","Ovary","Esophagus","Thymus","Lung","Liver","Pancreas","Stomach","Gastric","Intestine","Colon","Bladder"))
choose<-c()
for(i in 1:length(tissue)){
  tmp<-subset(data,group==tissue[i])
  tmp<-tmp[order(tmp[,3],decreasing=T),]
  if(nrow(tmp)>5){
  choose<-c(choose,tmp[1:80,1])  
  }else{
  choose<-c(choose,tmp[,1])    
  }
  choose
}
tmp2<-file1[match(choose,rownames(file1)),]
tmp2<-tmp2[,unlist(lapply(tissue,function(x) grep(x,colnames(tmp2))))]
write.table(tmp2,file="high.gsi.tissue.matrix.txt",sep='\t',quote=F,col.names=NA,row.names=T)

xx<-order(rowSums(tmp2),decreasing=T)[1:50]
yy<-which(unlist(apply(tmp2,1,function(x) max(x)<0.4)))
zz<-which(unlist(apply(tmp2,1,function(x) sum(x>0.5)>0.65*length(x))))

tmp2<-tmp2[-as.vector(c(xx,yy,zz)),]
tmp2[tmp2<0.4]<-0

library("grDevices")
library("gplots")
filename=paste("Figure-20-7-2-1-2-1",60,"pdf",sep=".")
pdf(filename)
col=colorRampPalette(c("yellow", "blue"))(20) 
heatmap.2(tmp2,col=col,trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.8,labRow=NA)
dev.off()


# tissues signatures
names(table(data[,2]))
lung.signature<-subset(data,GSI>0.5 & group=="Lung")              # 0.52 for bspp
colon.signature<-subset(data,GSI>0.55 & group=="Colon")
pancrease.signature<-subset(data,GSI>0.68 & group=="Pancreas")
nrow(lung.signature)
nrow(colon.signature)
nrow(pancrease.signature)

# prediction section
file1<-read.table("RRBS_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file1)
f2<-RawNARemove(file1,missratio=0.4)
file1<-impute.knn(data.matrix(f2))$data
# remove solid tissue
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])  # get sample name
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])        # get sample id
remove=c(samplename2[grep("6-T",samplename2)],samplename2[grep("PC-T",samplename2)],samplename2[grep("CTT-",samplename2)],samplename2[grep("7-T",samplename2)])
file1<-file1[,-match(remove,samplename2)]

samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
new<-read.table("/home/shg047/monod/phase2/saminfo.txt",sep="\t",as.is=T)
cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
matrix=file1
samplename2<-gsub("6-P","CC-P",samplename2)
samplename2<-gsub("7-P","LC-P",samplename2)
samplename2<-gsub("6-T","CC-T",samplename2)
samplename2<-gsub("7-T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
colnames(matrix)=samplename2

# random Forest
x<-t(matrix)
y<-as.factor(sapply(colnames(matrix),function(x) substr(x,1,2)))

fit<-randomForest(scale(x),y,importance=T)
top<-order(fit$importance[,6],decreasing=T)[1:500]
fit<-randomForest(scale(x)[,top],y,importance=T)
fit

bed1<-cor2bed(rownames(matrix))
bed2<-cor2bed(lung.signature[,1])
bed3<-cor2bed(colon.signature[,1])
bed4<-cor2bed(pancrease.signature[,1])

rlt.lung<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed2,opt.string="-wa -u")
rlt.colon<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed3,opt.string="-wa -u")
rlt.pancrease<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed4,opt.string="-wa -u")

cor.lung<-bed2cor(rlt.lung)
cor.colon<-bed2cor(rlt.colon)
cor.pancrease<-bed2cor(rlt.pancrease)

x<-y<-c()
z1<-z2<-z3<-z4<-z5<-c()
for(mhl in seq(0,0.3,by=0.01)){
  for(ratio in seq(0,0.2,by=0.01)){

# assess the prediction performance to cancer orign
data.prediction.lung<-file1[match(cor.lung,rownames(file1)),]                  # collect lung data
    
data.prediction.lung<-file1[match(cor.lung,rownames(file1)),c(grep("7-P",colnames(file1)))]                  # collect lung data
data.prediction.colon<-file1[match(cor.colon,rownames(file1)),c(grep("6-P",colnames(file1)))]                # collect colon data 
data.prediction.pancrease<-file1[match(cor.pancrease,rownames(file1)),c(grep("PC-P",colnames(file1)))]       # collect pancrease data
x1<-apply(data.prediction.lung,2,function(x) sum(x>mhl)/length(x))
x2<-apply(data.prediction.colon,2,function(x) sum(x>mhl)/length(x))
x3<-apply(data.prediction.pancrease,2,function(x) sum(x>mhl)/length(x))

# assess the prediction performance with Normal plasma as control
data.prediction.lung<-file1[match(cor.lung,rownames(file1)),grep("NC-P",colnames(file1))]
data.prediction.colon<-file1[match(cor.colon,rownames(file1)),grep("NC-P",colnames(file1))]
data.prediction.pancrease<-file1[match(cor.pancrease,rownames(file1)),grep("NC-P",colnames(file1))]
y1<-apply(data.prediction.lung,2,function(x) sum(x>mhl)/length(x))
y2<-apply(data.prediction.colon,2,function(x) sum(x>mhl)/length(x))
y3<-apply(data.prediction.pancrease,2,function(x) sum(x>mhl)/length(x))

print(c(x1,x2,x3,y1,y2,y3))
x<-c(x,mhl)
y<-c(y,ratio)
z1<-c(z1,sum(x1>ratio)/length(x1))                                        # accuracy for lung
z2<-c(z2,sum(x2>ratio)/length(x2))                                        # accuracy for colon
z3<-c(z3,sum(x3>ratio)/length(x3))                                        # accuracy for pancrease
z4<-c(z4,(sum(y1>ratio)+sum(y2>ratio)+sum(y3>ratio))/(3*length(y1)))      # accuracy for lung
z5<-c(z5,z1+z2+z3-z4)
  }
}
rlt1<-data.frame(x,y,z1,z2,z3,z4,z5)
write.table(rlt1,file="gsi.prediction.paramter.txt",quote=F,sep="\t")

### BSPP
file2<-read.table("/home/sguo/monod/phase2/150209_BSPP_mld_blocks_stringent_mhl_matrix.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file2)
samplename1=sapply(strsplit(colnames(file2),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
remove=c("6-T-3","6-T-4","7-T-2",paste("NC-P-",19:24,sep=""),"PC-P-10","6-P-6",paste("PC-P-",c(2,3,6,9),sep=""))
remove
file2<-file2[,-match(remove,samplename2)]
head(file2)
samplename1=sapply(strsplit(colnames(file2),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
new<-read.table("saminfo.txt",sep="\t",as.is=T)
cor1<-match(samplename2,new[,3])
lab1<-new[cor1,4]
groupname=lab1
matrix=file2
samplename2<-gsub("6-P","CC-P",samplename2)
samplename2<-gsub("7-P","LC-P",samplename2)
samplename2<-gsub("6-T","CC-T",samplename2)
samplename2<-gsub("7-T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
colnames(matrix)=samplename2
bed1<-cor2bed(rownames(matrix))
bed2<-cor2bed(lung.signature[,1])
bed3<-cor2bed(colon.signature[,1])
bed4<-cor2bed(pancrease.signature[,1])
rlt.lung<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed2,opt.string="-wa -u")
rlt.colon<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed3,opt.string="-wa -u")
rlt.pancrease<-Rbedtools(functionstring="intersectBed",bed1=bed1,bed2=bed4,opt.string="-wa -u")
cor.lung<-bed2cor(rlt.lung)
cor.colon<-bed2cor(rlt.colon)
cor.pancrease<-bed2cor(rlt.pancrease)
choose.prediction<-unique(c(cor.lung,cor.colon,cor.pancrease))
data.predition<-file2[match(choose.prediction,rownames(file2)),]
write.table(data.predition,file="BSPP.subset.WGBS.TSI.Predition.txt",sep="\t",col.names=NA,row.names=T,quote=F)  # send to desktop

prediction.data.plasma<-data.matrix(data.predition[,c(grep("6P",colnames(data.predition)),grep("7P",colnames(data.predition)),grep("NC",colnames(data.predition)))])
data.prediction.lung<-file1[match(cor.lung,rownames(file2)),c(grep("7P",colnames(file2)))]
data.prediction.colon<-file1[match(cor.colon,rownames(file2)),c(grep("6P",colnames(file2)))]
x1<-apply(data.prediction.lung,2,function(x) sum(x>0)/length(x))
x2<-apply(data.prediction.colon,2,function(x) sum(x>0/length(x)))

dim(data.prediction.lung)
dim(data.prediction.colon)

data.prediction.lung.normal.plasma<-file2[match(cor.lung,rownames(file2)),c(grep("NC",colnames(file2)))]
data.prediction.colon.normal.plasma<-file2[match(cor.colon,rownames(file2)),c(grep("NC",colnames(file2)))]
x1<-apply(data.prediction.lung.normal.plasma,2,function(x) sum(x>0)/length(x))
x2<-apply(data.prediction.colon.normal.plasma,2,function(x) sum(x>0)/length(x))

choose.prediction<-unique(c(cor.lung,cor.colon,cor.pancrease))
data.predition<-file1[match(choose.prediction,rownames(file1)),]
write.table(data.predition,file="RRBS.subset.WGBS.TSI.Predition.txt",sep="\t",col.names=NA,row.names=T,quote=F)  # send to desktop
prediction.data.plasma<-data.matrix(data.predition[,c(grep("6-P",colnames(data.predition)),grep("7-P",colnames(data.predition)),grep("PC-P",colnames(data.predition)))])
# take normal plasma as the control
prediction.data.plasma<-data.matrix(data.predition[,grep("NC-P",colnames(data.predition))])
# GSI estimation to prediction.plasma.data
colnames(prediction.data.plasma)<-c(rep("7-P",9),rep("6-P",10),rep("PC-P",5))
apply(prediction.data.plasma,2,function(x) sum(x>0)/length(x))

group=names(table(colnames(prediction.data.plasma)))
index=colnames(prediction.data.plasma)
gsi<-c()
gmaxgroup<-c()
for(i in 1:nrow(prediction.data.plasma)){
  gsit<-0
  gmax<-names(which.max(tapply(as.numeric(prediction.data.plasma[i,]),index,mean)))
  for(j in 1:length(group)){
    tmp<-(1-10^(mean(prediction.data.plasma[i,][which(index==group[j])]))/10^(mean(prediction.data.plasma[i,][which(index==gmax)])))/(length(group)-1)
    gsit<-gsit+tmp
  }
  gmaxgroup<-c(gmaxgroup,gmax)
  gsi<-c(gsi,gsit)
  print(c(gmax,gsit))
}
rlt=data.frame(region=rownames(prediction.data.plasma),group=gmaxgroup,GSI=gsi)
new.rlt<-rlt[order(rlt[,3],decreasing=T)[1:20],]

prediction.data.plasma.new<-prediction.data.plasma[match(new.rlt[,1],rownames(prediction.data.plasma)),]
library("grDevices")
library("gplots")
col=colorRampPalette(c("white", "red"))(20) 

pdf("Figure.RRBS.subset.WGBS.TSI.Predition.pdf")
prediction.data.plasma.new<-prediction.data.plasma.new[,order(colnames(prediction.data.plasma.new))]
heatmap.2(data.matrix(prediction.data.plasma.new),col=greenred(20),trace="none",density.info="none",Colv=F,Rowv=F,key=T,keysize=1,cexCol=0.7,labRow=NA)
dev.off()

x<-c(29,31,36,43,48,47,40,35,29,28)
names(x)<-c("-4*up","-3*up","-2*up","-1*up","Domain","-1*down","-2*down","-3*down","-4*down")
barplot(x,col="blue",ylim=c(0,50))


# Seq-Cap
file1<-read.table("/home/shg047/monod/oct/data/WGBS_SeqCap_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",row.names=1,as.is=T,check.names=F)
colnames(file1)
f2<-RawNARemove(file1,missratio=0.35)
file1<-impute.knn(data.matrix(f2))$data
# remove solid tissue
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])  # get sample name
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])        # get sample id
remove=c(samplename2[grep("6-T",samplename2)],samplename2[grep("PC-T",samplename2)],samplename2[grep("CTT-",samplename2)],samplename2[grep("7-T",samplename2)])
if(length(remove)>0){
  file1<-file1[,-match(remove,samplename2)]  
}else{
  file1<-file1
}
colnames(file1)
samplename1=sapply(strsplit(colnames(file1),"[.]"),function(x) unlist(x)[1])
samplename2=sapply(strsplit(samplename1,"_"),function(x) unlist(x)[1])
samplename2

samplename2<-gsub("6P","CC-P",samplename2)
samplename2<-gsub("7P","LC-P",samplename2)
samplename2<-gsub("6T","CC-T",samplename2)
samplename2<-gsub("7T","LC-T",samplename2)
samplename2<-gsub("frozen","Frozen",samplename2)
samplename2<-gsub("-100ng","",samplename2)
samplename2<-gsub("-5ng","",samplename2)
samplename2<-gsub("CTT","CC-T",samplename2)
samplename2<-gsub("PCP","PC-P",samplename2)
samplename2<-gsub("PCT","PC-T",samplename2)
samplename2<-gsub("NC-","NC-P-",samplename2)
samplename2
colnames(file1)=samplename2
colnames(file1)

x<-t(file1)
y<-as.factor(sapply(colnames(file1),function(x) substr(x,1,2)))
library("randomForest")
fit<-randomForest(scale(x),y,importance=T)
top<-order(fit$importance[,6],decreasing=T)[1:150]
fit<-randomForest(scale(x)[,top],y,importance=T)
fit
 
######## supplementary code
