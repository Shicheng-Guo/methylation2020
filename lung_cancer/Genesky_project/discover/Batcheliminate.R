BatchEvaluate<-function(expression_xls, sample_info_file, type='txt', write=T, covariates='all', par.prior=T, filter=F, skip=0, prior.plots=T){
	cat('Reading Sample Information Files\n')
	saminfo <- read.table(sample_info_file, header=T, sep='\t',comment.char='')
        if(sum(colnames(saminfo)=="Batch")!=1){return('ERROR: Sample Information File does not have a Batch column!')}
		
	cat('Reading Expression Data File\n')
	if(type=='csv'){
		dat <- read.csv(expression_xls,header=T,as.is=T)
                #print(dat[1:2,])
	#	dat <- dat[,trim.dat(dat)]  
                #print(colnames(dat))
		colnames(dat)=scan(expression_xls,what='character',nlines=1,sep=',',quiet=T)[1:ncol(dat)]
                #print(colnames(dat))
		}
         else{
		dat <- read.table(expression_xls,header=T,comment.char='',fill=T,sep='\t', as.is=T)
		dat <- dat[,trim.dat(dat)]
		colnames(dat)=scan(expression_xls,what='character',nlines=1,sep='\t',quiet=T)[1:ncol(dat)]
		}
	if (skip>0){
              geneinfo <- as.matrix(dat[,1:skip])
              dat <- dat[,-c(1:skip)]}
        else{geneinfo=NULL}
        #print(geneinfo[1:4])
        #print(dat[1:2,])
	
	if(filter){
		ngenes <- nrow(dat)
		col <- ncol(dat)/2
		present <- apply(dat, 1, filter.absent, filter)
		dat <- dat[present, -(2*(1:col))]
		if (skip>0){geneinfo <- geneinfo[present,]}
		cat('Filtered genes absent in more than',filter,'of samples. Genes remaining:',nrow(dat),'; Genes filtered:',ngenes-nrow(dat),'\n')
		}

	if(any(apply(dat,2,mode)!='numeric')){return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
	
	tmp <- match(colnames(dat),saminfo[,1])
	if(any(is.na(tmp))){return('ERROR: Sample Information File and Data Array Names are not the same!')}
	tmp1 <- match(saminfo[,1],colnames(dat))
	saminfo <- saminfo[tmp1[!is.na(tmp1)],]		

	if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
	design <- design.mat(saminfo)	

knn<-impute.knn(as.matrix(dat[,2:dim(dat)[2]]) ,k = 10, rowmax = 0.6, colmax = 0.6, maxp = 1500, rng.seed=362436069)
dat<-knn$data
dat<-normalize.quantiles(dat,copy=TRUE)
data<-data.frame(mydata[,1],data)
colnames(data)<-colnames(mydata)
write.table(data,file="GSE16559_GSE28094_TCGA_non-normalized_beta_data.knn.combat.txt",sep="\t",row.names=F,col.names=T,quote=F)
}


# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
	tmp <- strsplit(colnames(dat),'\\.')
	tr <- NULL
	for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
	tr
	}

# for GSE16559  
data<-read.table("GSE16559_non-normalized_data.txt",as.is=F,head=T,sep="\t")
i<-as.numeric();beta<-data[,1]
for(i in c(1:((dim(data)[2]-1)/2))){
   b<-data[,2*i]/(data[,2*i]+data[,2*i+1]+100)
   beta<-data.frame(beta,b)
}

sample<-colnames(data)
m<-colnames(data)[c(1,seq(2,length(colnames(data)),by=2))]
m<-gsub("\\.Signal\\.CY5","",m)

#a<-unlist(strsplit(sample,"\\."))
#b<-a[grep("X",a)]
#c<-unlist(strsplit(b,"X"))
#d<-((c[grep("\\d+",c,perl=T)]))
#e<-c(colnames(data)[1],d[seq(2,length(d),by=2)])
colnames(beta)<-m
write.table(beta,file="GSE16559_non-normalized_beta_data.txt",sep="\t",quote=F,row.names=F)
rm(list=ls())

#GSE28094_lung_non-normalized_data.txt
data<-read.table("GSE28094_lung_non-normalized_data.txt",as.is=F,head=T,sep="\t")

Pval<-data[,c(1,grep("Pval",colnames(data)))]
data<-data[,-grep("Pval",colnames(data))]
i<-as.numeric();beta<-data[,1]
for(i in c(1:((dim(data)[2]-1)/2))){
   b<-data[,2*i+1]/(data[,2*i]+data[,2*i+1]+100)
   beta<-data.frame(beta,b)
}

beta2<-array(as.matrix(beta[,seq(2,dim(beta)[2])]))
beta2[which(Pval[,seq(2,dim(Pval)[2])]>0.05)]<-NA
beta3<-matrix(beta2,dim(beta)[1])
beta4<-data.frame(beta[,1],beta3)

sample<-colnames(data)
m<-colnames(data)[c(1,seq(2,length(colnames(data)),by=2))]
m<-gsub("\\.CY3","",m)
colnames(beta4)<-m
write.table(beta4,file="GSE28094_lung_non-normalized_beta_data.txt",sep="\t",quote=F,row.names=F)
rm(list=ls())


#for GSE29505_unmethylated_and_methylated_signal_intensities.txt
data<-read.table("GSE29505_unmethylated_and_methylated_signal_intensities.txt",as.is=F,head=T,sep="\t")
i<-as.numeric();beta<-data[,1]
for(i in c(1:((dim(data)[2]-1)/2))){
   b<-data[,2*i+1]/(data[,2*i]+data[,2*i+1]+100)
   beta<-data.frame(beta,b)
}
sample<-colnames(data)
m<-colnames(data)[c(1,seq(2,length(colnames(data)),by=2))]
m<-gsub("\\.Unmethylated\\.Signal","",m)
colnames(beta)<-m
TargetID<-gsub("a\\_","",beta[,1])
TargetID<-gsub("\\-61","",TargetID)
beta[,1]<-TargetID
beta<-beta[,c(1,grep("Lung",colnames(beta)))]
write.table(colnames(beta),file="GSE29505.sample.txt",sep="\t",row.names=F,col.names=F)
write.table(beta,file="GSE29505_non-normalized_beta_data.txt",sep="\t",quote=F,row.names=F)

# For GPL
GPL<-read.table("goldengate.txt",as.is=F,head=T,sep="\t")
GPL27<-read.table("methy27k.anotation.txt",as.is=F,head=T,sep="\t")
match(as.factor(beta[,1]),GPL[,dim(GPL)[2]])

write.table(GPL,file="1.txt",sep="\t",quote=F,row.names=F)
write.table(as.factor(beta[,1]),file="2.txt",sep="\t",quote=F,row.names=F)

data4<-read.table("TCGA_Methylation.txt",sep="\t",head=T,as.is=F)
data5<-t(data4)
write.table(data5,file="TCGA_non-normalized_beta_data.txt",sep="\t",row.names=T,col.names=F)


#independent followed conde

data0<-read.table("GSE29505_non-normalized_beta_data.txt",sep="\t",head=T,as.is=F) #Feinburg  cgsite
data1<-read.table("GSE28094_lung_non-normalized_beta_data.txt",sep="\t",head=T,as.is=F) #GoldenGate  cgsite
data2<-read.table("TCGA_non-normalized_beta_data.txt",sep="\t",head=T,as.is=F) #methylation27k   cgsite
data3<-read.table("GSE16559_non-normalized_beta_data.txt",sep="\t",head=T,as.is=F) #GoldenGate  cgsite
x1<-match((data1[,1]),GPL1[,1])
x2<-match(GPL[,dim(GPL)[2]],GPL27[,1])
x3<-match((beta[,1]),GPL[,9])
m<-GPL27[x2,]
sum(!is.na(x2))
sum(!is.na(x3))
dim(data1)
dim(data2)
dim(data3)

#adjust every data to the first column is cgsite
x1<-match((data1[,1]),GPL1[,1])
data1[,dim(data1)[2]+1]<-data1[,1]
data1[,1]<-GPL1[x1,dim(GPL1)[2]]
x1<-match((data3[,1]),GPL1[,1])
data3[,dim(data3)[2]+1]<-data3[,1]
data3[,1]<-GPL1[x1,dim(GPL1)[2]]
colnames(data1)[1]
colnames(data2)[1]
colnames(data3)[1]
x1<-merge(data1,data2,by.x="TargetID",by.y="ID")
x2<-merge(x1,data3,by.x="TargetID",by.y="ID_REF")
write.table(x2,"GSE16559_GSE28094_TCGA_non-normalized_beta_data.txt",sep="\t",col.names=T,row.names=F)
match(data1[,1],data2[,1])
P121<-GPL27[x2[which(!is.na(x2))],]
dim(x1)
dim(x2)

write.table(P121,file="GoldenGate_Meth27K.share.probe.txt",sep="\t",row.names=F,col.names=F)
write.table(colnames(x2),file="sample.txt",sep="\t",row.names=F,col.names=F)
write.table(x2,file="lung_methyation_goldengate",sep="\t",row.names=F,col.names=T)
sample<-read.table("Sample.rdsf",sep="\t",head=T,as.is=F) #Feinburg  cgsite
data<-x2[,c(1,match(sample[,1],colnames(x2)))]
pdf("sample.pdf")
barplot(table(sample[,2]),col="green",main="Sample Number of each group",ylim=c(0,300))
dev.off()
write.table(data,file="GSE16559_GSE28094_TCGA_non-normalized_beta_data.txt",sep="\t",row.names=F,col.names=T)


para<-data.frame(data[,which(sample[,2]=="para")+1])
tumor<-data.frame(data[,which(sample[,2]=="tumor")+1])
normal<-data.frame(data[,which(sample[,2]=="normal")+1])
mydata<-data.frame(normal,para,tumor)

	library("affy")
	library("impute")
	library("preprocessCore")
	library("gplots")
	#remove columns which NAs over 30 percent
source("RM.R")
RM("GSE16559_GSE28094_TCGA_non-normalized_beta_data.txt")

mydata<-read.table("GSE16559_GSE28094_TCGA_non-normalized_beta_data.txt.rm",sep="\t",head=T,as.is=F) #Feinburg  cgsite
knn<-impute.knn(as.matrix(mydata[,2:dim(mydata)[2]]) ,k = 10, rowmax = 0.6, colmax = 0.6, maxp = 1500, rng.seed=362436069)
data<-knn$data
data<-normalize.quantiles(data,copy=TRUE)Sample_.txt
data<-data.frame(mydata[,1],data)
colnames(data)<-colnames(mydata)
write.table(data,file="GSE16559_GSE28094_TCGA_non-normalized_beta_data.knn.combat.txt",sep="\t",row.names=F,col.names=T,quote=F)

source("combat.R")
ComBat("GSE16559_GSE28094_TCGA_non-normalized_beta_data.knn.combat.txt","Sample_.txt",covariates="all",skip=1,filter=F,par.prior=T)

data<-read.table("Adjusted.GSE16559_GSE28094_TCGA_non-normalized_beta_data.knn.combat.txt",sep="\t",head=T,as.is=F,row.names=1) #Feinburg  cgsite
	fit1<-prcomp(t(data[,c(1:dim(data)[2])]),cor=TRUE,scale=TRUE)
	a1<-predict(fit1)

	variant<-read.table("Sample_.txt",as.is=F,head=T,sep="\t")

#PC1 Only
	# check sample type
        para<-which(variant[,2]=="para")
	tumor<-which(variant[,2]=="tumor")
	normal<-which(variant[,2]=="normal")
	pdf("PCA1.pdf")
	plot(a1[tumor,1],col="red",pch=1,xlab="Samples",ylab="Methylation PC1 of variance",xlim=c(-10,dim(a1)[1]+2),ylim=c(min(a1[,1]),max(a1[,1])+2))
	points(a1[para,1],col="blue",pch=2)
	points(a1[normal,1],col="green",pch=3)
        legend("topright",c("tumor","para","normal"),col=c("red","blue","green"),pch=c(1,2,3),cex=0.8)
	dev.off()
	# check batch effect[pca method]
        para<-which(variant[,3]=="1")
	normal<-which(variant[,3]=="3")
	tumor<-which(variant[,3]=="2")
	pdf("PCA2.pdf")
	plot(a1[tumor,1],col="red",pch=1,xlab="Samples",ylab="Methylation PC1 of variance",xlim=c(-10,dim(a1)[1]+2),ylim=c(min(a1[,1]),max(a1[,1])+2))
	points(a1[para,1],col="blue",pch=2)
	points(a1[normal,1],col="green",pch=3)
        legend("topright",c("Esteller","TCGA","Karl"),col=c("red","blue","green"),pch=c(1,2,3),cex=0.8)
	dev.off()
	# check batch effect[num method]
	factor<-




#PC1 & PC2
	# check sample type
        para<-which(variant[,2]=="para")
	tumor<-which(variant[,2]=="tumor")
	normal<-which(variant[,2]=="normal")
	pdf("PCA3.pdf")
	plot(a1[tumor,1],a1[tumor,2],col="red",pch=1,xlim=c(min(a1)-0.3*abs(min(a1)),1.2*max(a1)),ylim=c(min(a1)-0.3*abs(min(a1)),1.2*max(a1)),xlab="Methylation PC1 of variance",ylab="Methylation PC2 of variance")
	points(a1[para,1],a1[para,2],col="blue",pch=2)
	points(a1[normal,1],a1[normal,2],col="green",pch=3)
        legend("topright",c("tumor","para","normal"),col=c("red","blue","green"),pch=c(1,2,3),cex=0.8)
	dev.off()
	# check batch effect
        para<-which(variant[,3]=="1")
	tumor<-which(variant[,3]=="2")
	normal<-which(variant[,3]=="3")
	pdf("PCA4.pdf")
	
	plot(a1[tumor,1],a1[tumor,2],col="red",pch=1,xlim=c(min(a1)-0.3*abs(min(a1)),1.2*max(a1)),ylim=c(min(a1)-0.3*abs(min(a1)),1.2*max(a1)),xlab="Methylation PC1 of variance",ylab="Methylation PC2 of variance")
	points(a1[para,1],a1[para,2],col="blue",pch=2)
	points(a1[normal,1],a1[normal,2],col="green",pch=3)
        legend("topright",c("Esteller","TCGA","Karl"),col=c("red","blue","green"),pch=c(1,2,3),cex=0.8)
	dev.off()
	pdf("PCA2.me.1.pdf")
	plot(fit1,type="lines",main="AA")











