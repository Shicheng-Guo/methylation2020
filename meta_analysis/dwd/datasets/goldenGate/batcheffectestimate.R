
        ComBat("GSE16559_GSE28094_TCGA_non-normalized_beta_data.knn.combat.txt","Sample_.txt",covariates="all",skip=1,filter=F,par.prior=T)
        data<-read.table("GSE16559_GSE28094_TCGA_non-normalized_beta_data.knn.combat.txt",sep="\t",head=T,as.is=F) #Feinburg  cgsite
        data2<-read.table("Adjusted_GSE16559_GSE28094_TCGA_non-normalized_beta_data.knn.combat.txt",sep="\t",head=T,as.is=F) #Feinburg  cgsite

	fit1<-prcomp(t(data[,c(2:dim(data)[2])]),cor=TRUE,scale=TRUE)
	fit2<-prcomp(t(data2[,c(2:dim(data2)[2])]),cor=TRUE,scale=TRUE)
	a1<-predict(fit1)
	a2<-predict(fit2)
	variant<-read.table("Sample_.txt",as.is=F,head=T,sep="\t")



#PC1 Only
	# check sample type
	samplecolumn=2
	name<-names(table(variant[match(rownames(a1),variant[,1]),samplecolumn]))
	loc1<-which(variant[match(rownames(a1),variant[,1]),samplecolumn]==name[1])
	loc2<-which(variant[match(rownames(a1),variant[,1]),samplecolumn]==name[2])
	loc3<-which(variant[match(rownames(a1),variant[,1]),samplecolumn]==name[3])

	pdf("PCA1_1.pdf")	
	plot(a1[loc1,1],col="red",pch=1,xlab="Samples",ylab="Methylation PC1 of variance",xlim=c(-10,dim(a1)[1]+2),ylim=c(min(a1[,1]),max(a1[,1])+2))
	points(a1[loc2,1],col="blue",pch=2)
	points(a1[loc3,1],col="green",pch=3)
        legend("topright",name,col=1:length(name),pch=1:length(name),cex=0.8)
	dev.off()

	# check batch effect
	batchcolumn=3
	name<-names(table(variant[match(rownames(a1),variant[,1]),batchcolumn]))
	loc1<-which(variant[match(rownames(a1),variant[,1]),batchcolumn]==name[1])
	loc2<-which(variant[match(rownames(a1),variant[,1]),batchcolumn]==name[2])
	loc3<-which(variant[match(rownames(a1),variant[,1]),batchcolumn]==name[3])
	#loc4<-which(variant[match(rownames(a1),variant[,1]),batchcolumn]==name[4])
	#loc5<-which(variant[match(rownames(a1),variant[,1]),batchcolumn]==name[5])
	pdf("PCA1_2.pdf")
	plot(a1[loc1,1],col=1,pch=1,xlab="Samples",ylab="Methylation PC1 of variance",xlim=c(-10,dim(a1)[1]+2),ylim=c(min(a1[,1]),max(a1[,1])+2))
	points(a1[loc2,1],col=2,pch=2)
	points(a1[loc3,1],col=3,pch=3)
	#points(a1[loc4,1],col=4,pch=4)
	#points(a1[loc5,1],col=5,pch=5)
        legend("topright",name,col=1:length(name),pch=1:length(name),cex=0.8)
	dev.off()

#PC1 & PC2
	pdf("PCA2_2.pdf")
	par(mfrow=c(2,2))
	# check sample type
	name<-names(table(variant[match(rownames(a1),variant[,1]),samplecolumn]))
	loc1<-which(variant[match(rownames(a1),variant[,1]),samplecolumn]==name[1])
	loc2<-which(variant[match(rownames(a1),variant[,1]),samplecolumn]==name[2])
	loc3<-which(variant[match(rownames(a1),variant[,1]),samplecolumn]==name[3])

	plot(a1[loc1,1],a1[loc1,2],col=1,pch=1,xlab="Methylation PC1 of variance",ylab="Methylation PC2 of variance",xlim=c(-10,20),ylim=c(-10,10))
	points(a1[loc2,1],a1[loc2,2],col=2,pch=2)
	points(a1[loc3,1],a1[loc3,2],col=3,pch=3)
        legend("topright",name,col=1:length(name),pch=1:length(name),cex=0.8)

	loc21<-which(variant[match(rownames(a2),variant[,1]),samplecolumn]==name[1])
	loc22<-which(variant[match(rownames(a2),variant[,1]),samplecolumn]==name[2])
	loc23<-which(variant[match(rownames(a2),variant[,1]),samplecolumn]==name[3])

	plot(a2[loc21,1],a2[loc21,2],col=1,pch=1,xlab="Methylation PC1 of variance",ylab="Methylation PC2 of variance",xlim=c(-10,20),ylim=c(-10,10))
	points(a2[loc22,1],a2[loc22,2],col=2,pch=2)
	points(a2[loc23,1],a2[loc23,2],col=3,pch=3)
        legend("topright",name,col=1:length(name),pch=1:length(name),cex=0.8)


	# check batch effect
	name<-names(table(variant[match(rownames(a1),variant[,1]),3]))
	loc1<-which(variant[match(rownames(a1),variant[,1]),3]==name[1])
	loc2<-which(variant[match(rownames(a1),variant[,1]),3]==name[2])
	loc3<-which(variant[match(rownames(a1),variant[,1]),3]==name[3])
	#loc4<-which(variant[match(rownames(a1),variant[,1]),3]==name[4])
	#loc5<-which(variant[match(rownames(a1),variant[,1]),3]==name[5])

	loc21<-which(variant[match(rownames(a2),variant[,1]),3]==name[1])
	loc22<-which(variant[match(rownames(a2),variant[,1]),3]==name[2])
	loc23<-which(variant[match(rownames(a2),variant[,1]),3]==name[3])
	#loc24<-which(variant[match(rownames(a2),variant[,1]),3]==name[4])
	#loc25<-which(variant[match(rownames(a2),variant[,1]),3]==name[5])

	plot(a1[loc1,1],a1[loc1,2],col=1,pch=1,xlab="Methylation PC1 of variance",ylab="Methylation PC2 of variance",xlim=c(-10,20),ylim=c(-10,10))
	points(a1[loc2,1],a1[loc2,2],col=2,pch=2)
	points(a1[loc3,1],a1[loc3,2],col=3,pch=3)
	#points(a1[loc4,1],a1[loc4,2],col=4,pch=4)
	#points(a1[loc5,1],a1[loc5,2],col=5,pch=5)
        legend("topright",name,col=1:length(name),pch=1:length(name),cex=0.8)

	plot(a2[loc21,1],a2[loc21,2],col=1,pch=1,xlab="Methylation PC1 of variance",ylab="Methylation PC2 of variance",xlim=c(-10,20),ylim=c(-10,10))
	points(a2[loc22,1],a2[loc22,2],col=2,pch=2)
	points(a2[loc23,1],a2[loc23,2],col=3,pch=3)
	#points(a2[loc24,1],a2[loc24,2],col=4,pch=4)
	#points(a2[loc25,1],a2[loc25,2],col=5,pch=5)
        legend("topright",name,col=1:length(name),pch=1:length(name),cex=0.8)
	dev.off()

