
Spmarker<-function(mfile,pheno,tau=0.8,type='txt'){	
	library(gplots)
	library(impute)
#usage:
#source("Spmarker.R")
#Spmarker("mtcars.txt","saminfo.txt",tau=0.3)
#注意： 如果tau设的太大：会提示： `x' must be a numeric matrix
#saminfo.txt 中必须用Category作为标题,data的行为因子,列为个体

	cat('Reading Molecular Data File\n')
	dat <- read.table(mfile,header=T,comment.char='#',fill=T,sep='\t', as.is=T)
	cat('Reading Sample Information File\n')
	saminfo <- read.table(pheno, header=T, sep='\t',comment.char='#')
	if(sum(colnames(saminfo)=="Category")!=1){return('ERROR: Sample Information File does not have a Category column!')}

	dat<-impute.knn(as.matrix(dat[,2:ncol(dat)]))   
        Category<-as.numeric(saminfo$Category)
        data<-as.list(data.frame(dat$data,Category))
	aggre<-aggregate(data,by=list(Category),FUN=median)

	tau2<-as.numeric();
        for (i in 2:(ncol(aggre)-1)){
	tau1<-(nrow(aggre)-sum(aggre[,i])/max(aggre[,i]))/(nrow(aggre)-1)	
	tau2<-c(tau2,tau1)
	}

	color.map <- function(cancertype){ 
	cancertype
	}
	patientcolors <- as.character(unlist(lapply(saminfo$Category, color.map)))

	data<-as.data.frame(data)
	pdf("heatmap.pdf",width=17,height=17)  #increase it when fig is big
	heat<-heatmap.2(as.matrix(data[,which(tau2>tau)]),col=topo.colors(75), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors)
	dev.off()
	}



