RMNa<-function(expression_xls){
        data<-read.table(expression_xls,as.is=F,sep="\t",head=T,row.names=1)
	m1<-dim(data)[2];m2<-dim(data)[1];
	#remove columns which NAs over 30 percent
	a<-numeric();clen<-dim(data)[1];rlen<-dim(data)[2];i<-numeric();
	for (i in c(1:clen)){
		if (sum(is.na(data[,i]))>0.3*clen){  
			a<-c(a,i)
		}
	}
        if(sum(a)>0){
	data<-data[,-a];
	}
	print (paste("There are",dim(data)[2],"non-NA Cols at last;","Deltion ratio is" ,(m1-dim(data)[2])/m1,sep=" "))
	#remove rows which NAs over 30 percent
	b<-as.numeric()
	for (j in c(1:rlen)){
		if (sum(is.na(data[j,]))>0.3*rlen){  
		b<-c(b,j)
		}
	}
	if(sum(b)>0){
	data2<-data[-b,];   
	}
	ProbeID<-rownames(data2)
	print (paste("There are",dim(data2)[1],"non-NA Cols at last;","Deltion ratio is" ,(m2-dim(data2)[1])/m2,sep=" "))
	knn<-impute.knn(as.matrix(data2[,1:dim(data2)[2]]) ,k = 10, rowmax = 0.6, colmax = 0.6, maxp = 1500, rng.seed=362436069)
	dat<-knn$data
	dat<-normalize.quantiles(dat,copy=TRUE)
	colnames(dat)<-colnames(data2)
	ProbeID<-rownames(data2)	
	dat<-cbind(ProbeID,dat)
	outfile=paste(unlist(strsplit(expression_xls,"\\.txt"))[1],"narm.knn.quantitle.txt",sep=".");
	outfile2=paste(unlist(strsplit(expression_xls,"\\.txt"))[1],"t_narm.knn.quantitle.txt",sep=".");
	write.table(dat,file=outfile,sep="\t",quote=F,col.names=T,row.names=F)
	write.table(t(dat),file=outfile2,sep="\t",quote=F,col.names=T,row.names=F)
}

