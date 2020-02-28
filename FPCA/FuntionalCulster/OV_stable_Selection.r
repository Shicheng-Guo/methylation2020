#read the tables into the system
#the following code is saved under the new server:
#lma@HGCNT54:/home/lma/tcga_OV

chr <- read.table("FunPCAScoreGenome.txt.trans", header=T)

samp <- read.table("IndividualList.txt", header=T)
#extract the sample ID from the individual list
samp_ID <- data.frame(substr(samp[,1], 4,15))
colnames(samp_ID) <- c("ID")

#clinical information
clinical <- read.table("OV_Clinical.txt", header=T, sep='\t')
#the first and last columns that we need
cli <- cbind(data.frame(clinical[,1]), clinical[,9])
colnames(cli) <- c("ID","out")

#Merge them together
cli_merge <- merge(samp_ID, cli, all.x=T)

#extract the RNAseq data

chr_RNA <- chr[,which((cli_merge[,2] == "Resistant") |(cli_merge[,2]== "Sensitive"))]
RNA_resp <- cli_merge[which((cli_merge[,2] == "Resistant") |(cli_merge[,2]== "Sensitive")),]  

###############################################################
#split the sample into "Sensitive" and "Resistant"
chr_RNA_sen <- chr_RNA[,which(RNA_resp[,2] == "Sensitive")]
chr_RNA_Speci <- chr_RNA[,which(RNA_resp[,2] == "Resistant")] 

###############################################################
#draw 50% of the samples from the "Sen" and "Speci" respectively without replacement
save(chr_RNA_sen, chr_RNA_Speci, file="OV_RNA_split")

n.sen <- ncol(chr_RNA_sen)
n.speci <- ncol(chr_RNA_Speci)

#sample.sen <- chr_RNA_sen[,sample(n.sen,ceiling(n.sen/2))]
#sample.speci <- chr_RNA_Speci[,sample(n.speci,ceiling(n.speci/2))]
###############################################################
# Stable feature selection and classification algorithms
library(glmnet)
#assign the lambda values
	OV.out <- list()
	OV.sample <- list()
	for ( i in 1:100)
	{
	sample.sen <- chr_RNA_sen[,sample(n.sen,ceiling(n.sen/2))]
	sample.speci <- chr_RNA_Speci[,sample(n.speci,ceiling(n.speci*0.8))]
	#beacuse the sample of specificity is few

	#seed <- 123456
	new_y <- c(rep(1,ncol(sample.sen)),rep(0,ncol(sample.speci)))
	new_data <- cbind(sample.sen, sample.speci)
	
	OV.sample[[1]] <- colnames(new_data)
	
	######################################################
	x <- t(as.matrix(new_data))
	y <- matrix(new_y, ncol=1)

	#OV.las<-cv.glmnet(x, y,  family = "binomial", type="class")

	#OV.las<-glmnet(x, y,  family = "binomial", alpha=0.5)
	OV.las <- glmnet(x, y,  family = "binomial", lambda=seq(0.0001,0.5,0.005), alpha=0.5)
	#in this condition, if lambda > 0.35  then no variable will be selected
	
	#OV.las<-glmnet(x, y,  family = "binomial", alpha=0.5)

	#b <- as.matrix(coef(OV.las))
	#out <- rownames(b)[b[,100] != 0]
	OV.out[[i]] <- OV.las
	}

save(OV.out, OV.sample, file="OV.gla.NET.out")
######################################################
#summarize the output
#this is R code: 

# run.list <- list()
# beta.list <- list()
# lambda <- matrix(0,ncol=1, nrow=100)

# for (i in 1:100)
# {
	# b <- as.matrix(OV.out[[i]]$beta)	
	# for (j in 1:200)
	# {
		# out <- rownames(b)[b[,j] != 0]
		lambda[i,1] <- OV.out[[i]]$lambda[j]
		# beta.list[[j]] <- out
	# }
	# run.list[[i]] <- beta.list
# }

##################################################
#get the frequency of beta among 100 samplings

load("OV.gla.NET.out")

out.list <- list()

for (k in 1:100)
{
	beta.list <- list()
	
	for (i in 1:100)
	{
		#b <- as.matrix(OV.out[[i]]$beta)
		b <- as.matrix(OV.out[[i]]$beta)
		out <- rownames(b)[b[,k] != 0]
		#lambda[i,1] <- OV.out[[i]]$lambda[j]
		if (length(out) == 0)
		{
			beta.list[[i]] <- 0	
		}
		else
		{
			beta.list[[i]] <- out	
		}
	}

	beta.all <- matrix(0,1,1)
	rownames(beta.all) <- NULL
	
	for (ii in 1:100)
	{
			b <- as.matrix(beta.list[[ii]])		
			colnames(b) <- c("ID")
			
			if (b != 0)
			{
			#beta.all <- b
			
			beta <- merge(beta.all, b, all=T)
			beta.all <- beta
			}
			#beta.all <- beta
			
	}
	
##################################################
#calculate the frequency
	if (length(beta.all) != 0)
	{
	beta.freq <- matrix(0,nrow=nrow(beta.all), ncol=100)
	for (iii in 1:100)
	{
		b <- beta.list[[iii]]
	
		if(b !=0)
		{
		for (j in 1:nrow(beta.all))
		{
			if(beta.all[j,1] %in% b)
			{
				beta.freq[j,iii] <- 1
			}
		}
		
		}
	}
	
###################################################
	
	beta.sum <- matrix(0,nrow=nrow(beta.all), ncol=1)
	for (s in 1:nrow(beta.all))
	{
		beta.sum[s,1] <- sum(beta.freq[s,])
		rownames(beta.sum) <- beta.all[,1]
	}
	
	out.list[[k]] <- beta.sum
	}
}

save(out.list, OV.out, file="OV.stable.freq")
#####################################################
=====================================================
#####################################################
#organize the result:

select.gene <- NULL

for (i in 1:100)
{
	pp <- out.list[[i]]

	if (length(pp) !=0)
	{
	pp.rowname <- rownames(pp)
	tmp <-append(select.gene, pp.rowname)
	
	select.gene<- tmp
	}
}
##################
#and then get the unique gene list
select.mat <- unique(select.gene)

select.mat <- as.matrix(unique(select.gene))
result <- matrix(0,nrow=length(select.mat), ncol=1)

r <- cbind(select.mat,result)
colnames(r) <- c("ID", "freq")

#pp.o <- cbind(pp.out,as.matrix(rownames(pp.out)))
#colnames(pp.o) <- c("freq2", "ID")
#m.all <- merge(r, pp.o, all=T)
#colnames(pp.out) <- c("freq2")
#m.all <- merge(result,pp.out)


for (i in 1:100)
{
	pp <- out.list[[i]]

	if (length(pp) !=0)
	{
		#pp.out <- as.matrix(pp[rev(order(pp[,1])),])
		#pp.out.rowname <- rownames(pp.out)
		pp.out.rowname <- rownames(pp)
	
		#pp.mat <- cbind(pp.out.rowname, pp.out)
		pp.mat <- cbind(pp.out.rowname, pp)
		
		freq <- paste("order_", i, sep="")
		colnames(pp.mat) <- c("ID", freq)
	
		m.all <- merge(r, pp.mat, all=T)
		r <- m.all
	}
}

##########################################
#replace the NA values

rr <-r[,2]
for (i in 3:102)
{
	e <- as.numeric(as.character(r[,i]))
	e[is.na(e)] <- 0
	e1 <- as.matrix(e)
	mm <- cbind(rr, e1)
	rr <- mm
}
rownames(rr) <- r[,1]

rname <- rownames(rr)

new_name<- NULL;

for (i in 1:length(rname))
{
 new_name[i] <- strsplit(rname, "_", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[i]][1]
}

rownames(rr) <- new_name
	
kk <- rr[rev(order(rr[,101])),]

lambda <- OV.out[[1]]$lambda
final <- kk[1:300,2:101]

colnames(final) <- lambda

write.table(final, file="OV.RNASeq.FPCA.Stable.selection", col.names=T, row.names=T,sep="\t", quote=F)

