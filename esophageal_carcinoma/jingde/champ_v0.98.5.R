#############################champ
####written by Tiffany Morris
####updated 4 September 2013

load("probe_450K_VCs_af.RData")
load("probe_450K_features.RData")
load("probeInfoALL.Rd"); 

##################Upload, SVD and normalization#####################
champ.process<-function(rawdata=TRUE, directory = getwd(), resultsDir=paste(getwd(),"resultsChamp",sep="/"), methValue="B", filter = TRUE, detPcut=0.01, filterXY = TRUE, QCimages = TRUE, batchCorrect=T, runSVD = TRUE, studyInfo=F,norm = "BMIQ", runNorm=T, runDMR=T, runCNA=T,plotsBMIQ=T,plotCNA=T,bedFile=F, methProfile=F,controlProfile=F)
{
    batchDone=FALSE
    normDone=FALSE

	if(!file.exists(resultsDir))
	{
	
        (resultsDir)
		#message("Creating results directory. Results will be saved in ", resultsDir)
        #log<-print("Creating results directory. Results will be saved in ", resultsDir)
	
    }
	if(rawdata)
	{
		load=champ.load(directory=directory,methValue=methValue,resultsDir=resultsDir,QCimages=QCimages,filter=filter,detPcut=detPcut)
	}else
	{
		load=champ.read()
	}
	if(runNorm) #options are BMIQ, SWAN, PBC or NONE 
	{
		norm=champ.norm(beta=load$beta,rgSet=load$rgSet,pd=load$pd,mset=load$mset,methValue=methValue,plotsBMIQ=plotsBMIQ,QCimages=QCimages,norm=norm,resultsDir=resultsDir)
		normDone=T
	}
	
	if(runSVD)
	{
		if(normDone)
		{
			champ.SVD(beta=norm$beta,pd=load$pd,rgSet=load$rgSet,detP=load$detP, studyInfo=studyInfo,resultsDir=resultsDir)
		}else{
			champ.SVD(beta=load$beta,pd=load$pd,rgSet=load$rgSet,detP=load$detP, studyInfo=studyInfo,resultsDir=resultsDir)
		}	
	}
	if(batchCorrect)
	{
        if(normDone)
		{
            if(methValue=="B")
            {
                batchNorm=champ.runCombat(beta.c=norm$beta,pd=load$pd,logitTrans=TRUE)
            }else{
                batchNorm=champ.runCombat(beta.c=norm$beta,pd=load$pd,logitTrans=FALSE)
            }
            batchDone=TRUE
        }else{
            if(methValue=="B")
            {
                batchNorm=champ.runCombat(beta.c=load$beta,pd=load$pd,logitTrans=TRUE)
            }else{
                batchNorm=champ.runCombat(beta.c=load$beta,pd=load$pd,logitTrans=FALSE)
            }
            batchDone=TRUE
        }
	}
	
	if(runDMR)
	{
		if(batchDone==F & normDone)
		{
			dmr=champ.lasso(pd=load$pd,beta.norm=norm$beta,resultsDir=resultsDir,bedFile=bedFile,batchDone=batchDone)
		}else if(batchDone & normDone){
			dmr=champ.lasso(pd=load$pd,beta.norm=batchNorm$beta,resultsDir=resultsDir,bedFile=bedFile,batchDone=batchDone, normSave=norm$beta)
		}else if(batchDone==F & normDone==F)
        {
            dmr=champ.lasso(pd=load$pd,beta.norm=load$beta,resultsDir=resultsDir,bedFile=bedFile,batchDone=batchDone)
        }else if(batchDone & normDone == F)
        {
            dmr=champ.lasso(pd=load$pd,beta.norm=batchNorm$beta,resultsDir=resultsDir,bedFile=bedFile,batchDone=batchDone, normSave=load$beta)
        }
	}
    
	if(runCNA)
	{
		champ.CNA(intensity=load$intensity,pd=load$pd,batchCorrect=batchCorrect,resultsDir=resultsDir,plotCNA=plotCNA)
	}
	message("ChAMP complete.")
	
}

champ.read <- function(betaFile="beta.txt",sampleSheet="sampleSheet.txt",resultsDir)
{
	beta <-read.table(betaFile,row.names =T, sep = "\t")
	pd <- read.table(sampleSheet) #need to check
	if(file.exists(resultsDir))
	{
		message("The directory ",resultsDir," already exists. To avoid overwriting please rename in the arguments or delete this folder.")
		
	}else
	{ 
		dir.create(resultsDir)	
	}	
	return(list(beta=beta,pd=pd))
}

champ.load<-function(directory = getwd(),methValue="B",resultsDir=paste(getwd(),"resultsChamp",sep="/"),QCimages=TRUE,filter=TRUE,detPcut=0.01)
{
	if(require(minfi))
	{
	read.450k.sheet<-NA
	rm(read.450k.sheet)
	read.450k.exp<-NA
	rm(read.450k.exp)
	pData<-NA
	rm(pData)
	getGreen<-NA
	rm(getGreen)
	getRed<-NA
	rm(getRed)
	preprocessRaw<-NA
	rm(preprocessRaw)	
	getMeth<-NA
	rm(getMeth)
	getUnmeth<-NA
	rm(getUnmeth)		
	getBeta<-NA
	rm(getBeta)
	getM<-NA
	rm(getM)
	plotBetasByType<-NA
	rm(plotBetasByType)
	mdsPlot<-NA
	rm(mdsPlot)
    detectionP<-NA
    rm(detectionP)

	
	message("Loading data from ", directory)

	if(!file.exists(resultsDir))
	{
		dir.create(resultsDir)
		message("Creating results directory. Results will be saved in ", resultsDir)
	}			
		
    myDir= directory
    suppressWarnings(targets <- read.450k.sheet(myDir))
	rgSet <- read.450k.exp(base = myDir, targets = targets)
	sampleNames(rgSet)=rgSet[[1]]
	pd<-pData(rgSet)
	green=getGreen(rgSet)
	red=getRed(rgSet)
	mset <- preprocessRaw(rgSet)
    detP <- detectionP(rgSet)
	
    if(filter)
    {
        failed <- detP>detPcut
        numfail = colMeans(failed) # Fraction of failed positions per sample
        message("The fraction of failed positions per sample: ")
        print(numfail)
        fileName=paste(resultsDir,"/failedSample",".txt",sep="")
        write.table(numfail,fileName,row.names=T,col.names=paste("Sample_Name","Fraction_Failed_Samples",sep="\t"),quote=F,sep="\t")
        mset.f = mset[rowSums(detP <= detPcut) == ncol(detP),]
        message("By filtering with a detection p-value of ",detPcut," a total of ",dim(mset)[1]-dim(mset.f)[1]," probes have been removed from the analysis.")
        mset=mset.f
    }
        
    totalProbes=dim(mset)[1]
    intensity=getMeth(mset)+getUnmeth(mset)
    if(methValue=="B")
    {
        beta.raw = getBeta(mset, "Illumina")
    }else{beta.raw = getM(mset)}
    detP=detP[which(row.names(detP) %in% row.names(beta.raw)),]
    
	#cluster
	if(QCimages)
	{
		#save images...
        imageName1=paste(resultsDir,"raw_mdsPlot.pdf",sep="/")
        imageName2=paste(resultsDir,"raw_densityPlot.pdf",sep="/")
        imageName3=paste(resultsDir,"raw_SampleCluster.png",sep="/")
        main1=paste("Density plot of raw data (",totalProbes," probes)",sep="")
        main2=paste("All Samples before Normalization (",totalProbes, " probes)",sep="")
            
		pdf(imageName1,width=6,height=4)
		mdsPlot(mset,numPositions=1000,sampGroups=pd$Sample_Group,sampNames=pd$Sample_Name)
		dev.off()
				
		pdf(imageName2,width=6,height=4)
		densityPlot(rgSet,sampGroups=pd$Sample_Group,main=main1,xlab="Beta")
		dev.off()
		
		if(ncol(beta.raw) < 60)
		{	
			png(imageName3)
			betar_d<-dist(t(beta.raw))
			plot(hclust(betar_d),main=main2,cex=0.8)
			dev.off()
		}else{
			message("Cluster image is not saved when the number of samples exceeds 60.")
			}
	}
	return(list(mset=mset,rgSet=rgSet,pd=pd,intensity=intensity,beta=beta.raw,detP=detP));
}}

champ.colname<-function(beta,pd)
{
	c=colnames(beta)
	for( i in 1:length(c)){c[i]= pd$Sample_Name[i]}
	colnames(beta)=c
	return(beta)
}

champ.norm<-function(beta=load$beta,rgSet=load$rgSet,pd=load$pd,mset=load$mset,sampleSheet="sampleSheet.txt",resultsDir=paste(getwd(),"resultsChamp",sep="/"),methValue="B",fromIDAT=T, norm="BMIQ", fromFile = F, betaFile, filterXY=F, QCimages=T,plotsBMIQ=T)
{
	detectionP<-NA
	rm(detectionP)
	preprocessSWAN<-NA
	rm(preprocessSWAN)
	getBeta<-NA
	rm(getBeta)	
	getM<-NA
	rm(getM)
	cwd=getwd()
	#probe.features<-NULL
	#probeInfoALL.lv<-NULL			

	if(require(minfi)){	
	message("Normalizing data with ",norm)
	if(fromIDAT==F)
	{
		if(norm=="SWAN")
		{
			message("Swan normalization is not available from this data. You need to use champ.load() with raw IDAT files.")
			return()
			
		}else if(fromFile==F)
		{
			beta.p=betaFile
		}
		else
		{
			#check if genome studio file
			beta.p=champ.read(betaFile,sampleSheet)
		}
	}else
	{
        if(norm == "SWAN")
		{
			mset <-preprocessSWAN(rgSet, mset)
			if(methValue=="B")
			{	
				beta.p = getBeta(mset, "Illumina")
			}else{beta.p = getM(mset)}
		}else{
        beta.p=beta
        }
	}	
	if(norm=="BMIQ" | norm == "PBC")
	{
		### create design.v file
		match(rownames(beta.p),probeInfoALL.lv[[5]]) -> map.idx;
		mapto <- function(tmp.v){ return(tmp.v[map.idx]);}
 		probeInfo.lv  <- lapply(probeInfoALL.lv,mapto);
		design.v <- probeInfo.lv[[2]];

		if(norm == "BMIQ")
 		{
			newDir=paste(resultsDir,"Normalization",sep="/")
 			if(plotsBMIQ){if(!file.exists(resultsDir)){dir.create(resultsDir)}
 			if(!file.exists(newDir))
			{
				dir.create(newDir)
			}}
 			design.v<-as.numeric(design.v)
 			bmiq=beta.p
			hf.v <- vector();

			if(plotsBMIQ){setwd(newDir)}
 			for(s in 1:ncol(beta.p))
 			{
				sID=colnames(beta.p)[s]
				beta.v <- beta.p[,s];


				bmiq.o <- BMIQ(beta.v,design.v,doH=TRUE,nL=3,nfit=5000,niter=10,plots=plotsBMIQ,sampleID=sID);
				bmiq[,s] <- bmiq.o$nbeta;
				hf.v[s] <- bmiq.o$hf;
				
				print(sID);
			}
			setwd(cwd)
			beta.p=bmiq
		}else if(norm == "PBC")
		{
			beta.p=DoPBC(beta.p,design.v)
		}
	}else if(norm == "NONE")
	{
		
	}else{
		if(norm!="SWAN"){
		message("You have not selected a valid normalization method. No normalization will be preformed.")}
		}
	#else if(norm == "SubQ")
	#{
		#beta.norm=DoSubQ(beta.raw)
		#add normalization to colname _subQ

	#}
	if(filterXY)
	{
		autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
		beta.p=beta.p[row.names(beta.p) %in% autosomes$ID, ]	
	}
    totalProbes=dim(beta.p)[1]
	if(QCimages)
	{
		#save images...
		imageName=paste(resultsDir,"norm_mdsPlot.pdf",sep="/")
		pdf(imageName,width=6,height=4)
		mdsPlot(beta.p,numPositions=1000,sampGroups=pd$Sample_Group,sampNames=pd$Sample_Names)
		dev.off()
		
		imageName=paste(resultsDir,"norm_densityPlot.pdf",sep="/")
		pdf(imageName,width=6,height=4)
		densityPlot(beta.p,sampGroups=pd$Sample_Group,main="Density plot of normalized data",xlab=methValue)
		dev.off()
		
        if(ncol(beta.p) < 60)
        {
            #cluster
            imageName=paste(resultsDir,"norm_SampleCluster.png",sep="/")
            png(imageName)
            betar_d<-dist(t(beta.p))
            plot(hclust(betar_d),main="All Samples after Normalization",cex=0.8)
            dev.off()
    }

	}
	message("The analysis will proceed with ", dim(beta.p)[1], " probes and ",dim(beta.p)[2], " samples.")
	return(list(beta=beta.p))
}}

###########SVD##############
champ.SVD<-function(beta=norm$beta, rgSet=load$rgSet,detP=load$detP,pd=load$pd, loadFile=F, betaFile="beta.txt", sampleSheet="sampleSheet.txt", methProfile=F, methFile="MethylationProbeProfile.txt", controlProfile=F, controlFile="ControlProbeProfile.txt", studyInfo=F,studyInfoFile="studyInfo.txt",resultsDir=paste(getwd(),"resultsChamp",sep="/"))
{
    impute.knn<-NA
	rm(impute.knn)
    
	if(require(impute)){if(require(marray)){if(require(minfi)){
	message("Performing SVD")
	if(!is.null(rgSet))
	{
		if(is.null(detP)){detP <- detectionP(rgSet)}
		if(is.null(beta))
		{
			mset <- preprocessRaw(rgSet)
			beta=getBeta(mset,"Illumina")
			
		}
	}else {
		if(is.null(detP))
		{
			message("You need a matrix of detection p-values to run SVD. Please rerun with parameter detP")
			return()
		}
	}
	
	if(loadFile)
	{
		beta.p=champ.read(betaFile,sampleSheet)
		beta.m=beta.p$beta
		pd=beta.p$pd
	}else
	{
		beta.m<-beta		
	}
	if(methProfile)
	{
		if(!file.exists(methFile)){message("You don't have a MethylationProbeProfile.txt file in this directory"); return()}
		tmp.df<- read.table(methFile,header=TRUE,sep="\t",fill=T);
		tmpTarget <- as.vector(tmp.df$TargetID);
		tmp.df<- tmp.df[,c(grep("X(.*)\\.",names(tmp.df)))]

		tmpB.idx <- grep("AVG_Beta",colnames(tmp.df));
		tmpI.idx <- grep("Intensity",colnames(tmp.df));
		tmpU.idx <- grep("Signal_A",colnames(tmp.df));
		tmpM.idx <- grep("Signal_B",colnames(tmp.df));
		tmpPV.idx <- grep("Detection",colnames(tmp.df));

		data.l <- list(B=as.matrix(tmp.df[,tmpB.idx]),I=as.matrix(tmp.df[,tmpI.idx]),U=as.matrix(tmp.df[,tmpU.idx]),M=as.matrix(tmp.df[,tmpM.idx]),pv=as.matrix(tmp.df[,tmpPV.idx]));
		for(i in 1:length(data.l))
		{
  			rownames(data.l[[i]]) <- tmpTarget;
		}
		rm(tmp.df,tmpB.idx,tmpI.idx,tmpU.idx,tmpM.idx,tmpPV.idx,tmpTarget);

		for(i in 1:length(data.l))
		{
  			tmp.v <- gsub(".AVG_Beta","",colnames(data.l[[i]]));
  			tmp.v <- gsub(".Detection.Pval","",tmp.v);
 			tmp.v <- gsub(".Intensity","",tmp.v);
  			tmp.v <- gsub(".Signal_A","",tmp.v);
  			tmp.v <- gsub(".Signal_B","",tmp.v);
  			colnames(data.l[[i]]) <- tmp.v;
		}
		beta.m <- data.l$B
		rm(tmp.v)
	}

	if(!is.null(rgSet))
	{
		bc1=getControlAddress(rgSet,controlType=c("BISULFITE CONVERSION I"))
		bc2=getControlAddress(rgSet,controlType=c("BISULFITE CONVERSION II"))
		ext=getControlAddress(rgSet,controlType=c("EXTENSION"))
		tr=getControlAddress(rgSet,controlType=c("TARGET REMOVAL"))
		hyb=getControlAddress(rgSet,controlType=c("HYBRIDIZATION"))
			
		#columns are samples
		CPP=rbind(getGreen(rgSet)[bc1[1:3],],getRed(rgSet)[bc1[7:9],],getRed(rgSet)[bc2[1:4],],getGreen(rgSet)[tr[1:2],],getGreen(rgSet)[hyb[1:3],],getRed(rgSet)[ext[1:2],],getGreen(rgSet)[ext[3:4],])
		controlNames <- c("BSC-I C1 Grn","BSC-I C2 Grn","BSC-I C3 Grn","BSC-I C4 Red","BSC-I C5 Red","BSC-I C6 Red","BSC-II C1 Red","BSC-II C2 Red","BSC-II C3 Red","BSC-II C4 Red","Target Removal 1 Grn","Target Removal 2 Grn","Hyb (Low) Grn","Hyb (Medium) 		Grn","Hyb (High) Grn","Extension (A) Red","Extension (T) Red","Extension (C) Grn","Extension (G) Grn");
		rownames(CPP)=controlNames
			
		#log2
		dataC2.m <- matrix(nrow=length(controlNames),ncol=ncol(beta.m));
		colnames(dataC2.m) <- colnames(beta.m);
		rownames(dataC2.m) <- controlNames;
		for(r in 1:nrow(dataC2.m))
		{
  			dataC2.m[r,] <- log2(as.numeric(CPP[r,]));
		}

	}else{	
		if(controlProfile)
		{	
			
		ControlProbeProfiles<-read.delim(controlFile)
		ControlProbeProfiles_s <- ControlProbeProfiles[,grep("Signal",colnames(ControlProbeProfiles))];
		ControlProbeProfiles_t <- ControlProbeProfiles[,grep("ID",colnames(ControlProbeProfiles))];
		ControlProbeProfiles <- cbind(ControlProbeProfiles_s,ControlProbeProfiles_t);

		dataC.m <- ControlProbeProfiles;

		dataC.m<- dataC.m[,c(grep("X(.*)\\.",names(dataC.m)))]
		dataC.l <- list();
		dataC.l[[1]] <- dataC.m[,grep("Signal_Grn",colnames(dataC.m))];
		dataC.l[[2]] <- dataC.m[,grep("Signal_Red",colnames(dataC.m))];
		names(dataC.l) <- c("Grn","Red");
		colnames(dataC.l$Grn) <- gsub("X","",gsub("_","",gsub(".Signal_Grn","",colnames(dataC.l$Grn))));
		colnames(dataC.l$Red) <- gsub("X","",gsub("_","",gsub(".Signal_Red","",colnames(dataC.l$Red))));

		rownames(dataC.l$Grn) <- ControlProbeProfiles$ProbeID
		rownames(dataC.l$Red) <- ControlProbeProfiles$ProbeID
		
		#if(studyInfo)
		#{
		#	clin.df <- read.csv(studyInfoFile,header=TRUE,fill=T,sep="\t");
		#	pos1 <- pd$Sentrix_ID
		#	pos2 <- pd$Sentrix_Position
		#	sampleID.v <- paste(as.vector(pd$Sentrix_ID),"_",as.vector(pd$Sentrix_Position),sep="");
		#	clinID.v <- as.vector(pd$Sample_Names);

		#	match(colnames(data.l[[1]]),clinID.v) -> map.idx;
		#	clin.m <- as.matrix(clin.df[map.idx,2:ncol(clin.df)]);
		#	rownames(clin.m) <- as.vector(clin.df[map.idx,1]); 
        #ß}

		### build probe control data matrix
		seltypeC.v <- c("BSC-I C1 Grn","BSC-I C2 Grn","BSC-I C3 Grn","BSC-I C4 Red","BSC-I C5 Red","BSC-I C6 Red","BSC-II C1 Red","BSC-II C2 Red","BSC-II C3 Red","BSC-II C4 Red","Target Removal 1 Grn","Target Removal 2 Grn","Hyb (Low) Grn","Hyb (Medium) Grn","Hyb (High) Grn","Extension (A) Red","Extension (T) Red","Extension (C) Grn","Extension (G) Grn");
		seltypeC.idx <- c(1:3,7:9,13:16,833:834,21:23,17:20); #this chooses the specific control columns...
		grn.idx <- grep("Grn",seltypeC.v);
		red.idx <- grep("Red",seltypeC.v);

		## select specific controls
		dataC.m <- matrix(nrow=length(seltypeC.v),ncol=ncol(beta.m));
		colnames(dataC.m) <- colnames(beta.m);
		rownames(dataC.m) <- seltypeC.v;

		for(i in 1:ncol(dataC.m))
		{
			dataC.m[grn.idx,][,i] <- dataC.l$Grn[seltypeC.idx[grn.idx],][,i];
			dataC.m[red.idx,][,i] <- dataC.l$Red[seltypeC.idx[red.idx],][,i];
		}

		#log2
		dataC2.m <- matrix(nrow=length(seltypeC.v),ncol=ncol(beta.m));
		colnames(dataC2.m) <- colnames(beta.m);
		rownames(dataC2.m) <- seltypeC.v;
		for(r in 1:nrow(dataC2.m))
		{
  			dataC2.m[r,] <- log2(as.numeric(dataC.m[r,]));
		}
		dataC2.m=t(dataC2.m)
		}
	}

	PhenoTypes.lv <- list();

	################Customise Phenotype Data########################
	pd=pd[order(pd$Sample_Name),]
	well=unique(pd$Sample_Well)
	plates=unique(pd$Sample_Plate)
	groups=unique(pd$Sample_Group)
	pool=unique(pd$Pool_ID)
	chips=unique(pd$Slide)
	arrays=unique(pd$Array)


	allP=list(c(well,plates,groups,pool,chips,array))

	p=1
	###Well
	if(length(well)>1)
	{
		tmp.v <- as.vector(pd$Sample_Well);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v)); 
		for(i in 2:length(well))
		{
			PhenoTypes.lv[[p]][grep(well[i],tmp.v)] <- i; 
		}
		names(PhenoTypes.lv)[p]<-"Sample_Well"
		p=p+1
	}
	
	##Plate
	if(length(plates)>1)
	{
		tmp.v <- as.vector(pd$Sample_Plate);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v)) #1
		for(i in 2:length(plates))
		{
			PhenoTypes.lv[[p]][grep(plates[i],tmp.v)] <- 2
		}
		names(PhenoTypes.lv)[p]<-"Sample_Plate"
		p=p+1
	}

	##Group
	if(length(groups)>1)
	{
		tmp.v <- as.vector(pd$Sample_Group);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v));
		for(i in 2:length(groups))
		{
			PhenoTypes.lv[[p]][grep(groups[i],tmp.v)] <- i; #R
		}
		names(PhenoTypes.lv)[p]<-"Sample_Group"
		p=p+1
	}
	
	##Pool_ID
	if(length(pool)>1)
	{
		tmp.v <- as.vector(pd$Pool_ID);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v));
		for(i in 2:length(pool))
		{
			PhenoTypes.lv[[p]][grep(pool[i],tmp.v)] <- i;
		}
		names(PhenoTypes.lv)[p]<-"Pool_ID"
		p=p+1
	}

	### Slide
	if(length(chips)>1)
	{
		tmp.v <- as.vector(pd$Slide);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v));
		#check for just 1
		for(i in 2:length(chips))
		{
			PhenoTypes.lv[[p]][grep(as.numeric(chips[i]),tmp.v)] <- i; 
		}
		names(PhenoTypes.lv)[p]<-"Slide"
		p=p+1
	}

	### Sentrix_Position
	if(length(arrays)>1)
	{
		tmp.v <- as.vector(pd$Array);
		PhenoTypes.lv[[p]] <- rep(1,length(tmp.v)); ## R01C01
		for(i in 2:length(arrays))
		{
			PhenoTypes.lv[[p]][grep(arrays[i],tmp.v)] <- i; ## R01C02
		}
		names(PhenoTypes.lv)[p]<-"Array"
		p=p+1
	}
	########add in studyInfo

	if(studyInfo)
	{
		tmp.info<- read.table(studyInfoFile,header=TRUE,sep="\t",fill=T,stringsAsFactors=FALSE,as.is=T);
		tmp.info=tmp.info[which(tmp.info$Sample_Name %in% pd$Sample_Name),]
		if(!is.null(tmp.info))
		{
			tmp.info=tmp.info[order(tmp.info$Sample_Name),]
			tmp.info=as.matrix(tmp.info)
			#match tmp.info[1] with pd=$Sample_Names
			for(n in 2:ncol(tmp.info))
			{	
				p=p+1
				tmp.v <- as.vector(tmp.info[,n]);
				if(length(unique(tmp.v))>1)
				{
					PhenoTypes.lv[[p]] <- rep(1,length(tmp.v))
					phenos=unique(tmp.info[,n])
					for(i in 2:length(phenos))
					{
						PhenoTypes.lv[[p]][grep(phenos[i],tmp.v)]<-i
					}
					names(PhenoTypes.lv)[p]<-colnames(tmp.info)[n]
				}else(p=p-1)
			}
		}else{ message("Your sample info file doesn't have the same samples as the sample sheet")}
	}
	p=p-1
	factor.log <- c(rep(TRUE,p))
	if(!is.null(rgSet) | controlProfile==T)
	{
		factor.log <- c(factor.log,rep(FALSE,14))
	}

	##DoQC.R
	common.idx <- 1:nrow(beta.m);
	selCL.idx <- 1:ncol(beta.m);

	### find coverage per sample
    ndet.v <- vector();
	det.li <- list();

	#detection p-value
    for(s in 1:ncol(beta.m))
	{
		det.li[[s]] <- which(detP[,s]<0.05);
	  	ndet.v[s] <- length(det.li[[s]]);
	}
	covS.v <- ndet.v/nrow(beta.m);
	
    ### global coverage
	for(s in selCL.idx)
	{
		common.idx <- intersect(common.idx,det.li[[s]]);
	}

	### impute remaining missing values

    capture.output(irawdataCL.m <- impute.knn(beta.m[common.idx,selCL.idx],k=5)$data)
	
    ### DoSVD.R
    tmp.m<- irawdataCL.m-rowMeans(irawdataCL.m)
	rmt.o <- EstDimRMTv2(tmp.m);
	svd.o <- svd(tmp.m);
    if(rmt.o$dim >6)
    {
        topPCA <- 6;
    }else{topPCA <- rmt.o$dim}
        
	#Check coloumn idx & change to appropriate coloumns from dataC2.m(below includes control profile, sentrix and clin.info)
    selcat.idx <- c(1:length(PhenoTypes.lv));
	if(!is.null(rgSet) | controlProfile)
	{
		svdPV.m <- matrix(nrow=topPCA,ncol=length(selcat.idx)+nrow(dataC2.m));
		colnames(svdPV.m) <- c(names(PhenoTypes.lv)[selcat.idx],rownames(dataC2.m));
		print(cbind(factor.log[1:length(selcat.idx)],names(PhenoTypes.lv)[selcat.idx]));
	}else
	{
		svdPV.m <- matrix(nrow=topPCA,ncol=length(selcat.idx));
		colnames(svdPV.m) <- c(names(PhenoTypes.lv)[selcat.idx]);
		print(cbind(factor.log[1:length(selcat.idx)],names(PhenoTypes.lv)[selcat.idx]));
	}
	for(c in 1:topPCA)
	{
		tmp.v <- svd.o$v[,c];
	  	for(f in 1:length(selcat.idx))
	  	{
    		if(factor.log[f])
    		{
    	  		svdPV.m[c,f] <- kruskal.test(tmp.v ~ as.factor(PhenoTypes.lv[[selcat.idx[f]]]))$p.value;
    		}else 
    		{
      			svdPV.m[c,f] <- summary(lm(tmp.v ~ PhenoTypes.lv[[f]]))$coeff[2,4];
    		}
  		}
  		if(!is.null(rgSet) | controlProfile)
  		{
  			for(f in (length(selcat.idx)+1):ncol(svdPV.m))
  			{
    			svdPV.m[c,f] <- summary(lm(tmp.v ~ dataC2.m[f-length(selcat.idx),selCL.idx]))$coeff[2,4];
  			}
  		}
  		print(c);
	}
	### image heatmap
	myPalette <- c("darkred","red","orange","pink","white");
	breaks.v <- c(-200,-10,-5,-2,log10(0.05),0);
	imageName=paste(resultsDir,"SVDsummary.pdf",sep="/")
	pdf(imageName,width=8,height=8);
	par(mar=c(5,15,2,1),xpd=TRUE);
	image(x=1:nrow(svdPV.m),y=1:ncol(svdPV.m),z=log10(svdPV.m),col=myPalette,breaks=breaks.v,xlab="",ylab="",axes=FALSE, main= "Singular Value Decomposition Analysis (SVD)");
	axis(1,at=1:nrow(svdPV.m),labels=paste("PC-",1:nrow(svdPV.m),sep=""),las=2);
	suppressWarnings(axis(2,at=1:ncol(svdPV.m),labels=colnames(svdPV.m),las=2));
    legend(x=-3, y=5, legend=c(expression("p < 1x"~10^{-10}), expression("p < 1x"~10^{-5}),"p < 0.01", "p < 0.05", "p > 0.05"), fill=c("darkred","red","orange","pink","white"));    
	dev.off();

	#imageName=paste(resultsDir,"ScatterPC12.pdf",sep="/")
	#pdf(imageName,width=10,height=10);
	#plot(svd.o$v[,1],svd.o$v[,2],col=PhenoTypes.lv$Sample_Group,pch=c(15,15,15)[PhenoTypes.lv$Group],xlab="PC1",ylab="PC2",main="Plot of Priciple Components 1 & 2 showing groups");
	#legend('topright',c(groups[1],groups[2]),fill=c("black","red"),bty="n");
	#dev.off();
    }}}}

####################Batch Correct################
champ.runCombat <- function(beta.c=norm$beta,pd=load$pd,logitTrans=TRUE,batch=pd$Slide)
{	
	if(require(sva)){		
	message("Beginning batchCorrect")
 	message("Preparing files for ComBat")
 	mod = model.matrix(~as.factor(Sample_Group), data=pd)
 	
    if(length(unique(batch))<2)
	{
		message("You must have at least two different slides for Combat to correct for batch effects")
		return(list(beta=beta.c,pd=pd))
	}else{

        if(logitTrans)
        {
            message("Your data will be logit transformed, batch corrected and inverse logit transformed")
            log=logit2(beta.c)
            combat=champ.ComBat(dat=log,batch=batch,mod=mod,par.prior=TRUE)
            if(!is.null(combat)){combat=ilogit2(combat)}
        }else{
            message("Your data will be batch corrected")
            combat=champ.ComBat(dat=beta.c,batch=batch,mod=mod,par.prior=TRUE)

        }
        if(is.null(combat))
        {
            return(list(beta=beta.c,pd=pd))
        }else{
            return(list(beta=combat,pd=pd))
        }
	}
}}

###########MVP##############
champ.MVP<-function(beta.norm = norm$beta, pd=load$pd, adjPVal=0.05, adjust.method="BH", compare.group=c("C","T"),resultsDir=paste(getwd(),"resultsChamp",sep="/"),bedFile=TRUE)
{
	#put them in order based on data
	#limma
	#make a bedfile and then call DMR function.
	
	if(require(limma)){
	makeContrasts<-NA
	rm(makeContrasts)
	lmFit<-NA
	rm(lmFit)
	contrasts.fit<-NA
	rm(contrasts.fit)
	eBayes<-NA
	rm(eBayes)
	topTable<-NA
	rm(topTable)	
	
    groupLabel=unique(pd$Sample_Group)
	if(length(groupLabel)>2)
	{message("Your dataset has more than two groups. Analysis will compare the first two")}
	
	controls=pd[which(pd$Sample_Group==groupLabel[1]),]
	test=pd[which(pd$Sample_Group==groupLabel[2]),]
	all=c(controls$Sample_Name,test$Sample_Name)
	data=matrix(NA,length(row.names(beta.norm)),length(all))
	row.names(data)=row.names(beta.norm)
	colnames(data)=all
	
	#quicker way to do this?
	if(length(beta.norm)==0)
	{
		message("Your dataset is empty and MVP list cannot be produced")
		return()
	}
	for(i in 1:length(all))
	{
		done=F
		j=1
		while(done==F)
		{
			if(colnames(data)[i]==colnames(beta.norm)[j])
			{
				data[,i]=beta.norm[,j] 
				done=T
			}
			else
			{
				if(j<ncol(beta.norm))
				{
					j=j+1
				}
				else
				{
					print("There is an error in this dataset")
				}
			}
		}
		
	}

	numsamples=length(controls$Sample_Name)
	numtest=length(test$Sample_Name)
	design <- model.matrix(~0 + factor(c(rep("C",numsamples), rep("T",numtest))))
	print(paste("contrast",groupLabel[1],groupLabel[2],sep=" "))
	colnames(design) <- c("C", "T")
	contrast.matrix <- makeContrasts(T-C, levels=design)

	fit <- lmFit(data, design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	ok=TRUE
	tryCatch(fit3 <- eBayes(fit2),
      warning=function(w) 
      {
      	cat("No sample variance.\n")
      	ok <- FALSE
      })
      if(ok)
      {
          results<- topTable(fit3, coef=1, number=dim(data)[1], adjust.method=adjust.method,p.value=adjPVal)
          message("You have found ", dim(results)[1], " significant MVPs with a ",adjust.method," adjusted P-value below ", adjPVal)
          if(dim(results)[1]==0)
          {
              message("No bedfile will be generated for tophits but a full MVP list with all p-values is being saved")
          }else{
          resList=data[which(row.names(data) %in% results$ID),]
          resList=data.frame(resList)
	
          #join with annotation
          resList_anno <- data.frame(resList, probe.features[match(rownames(resList), rownames(probe.features)),])
          
          if(bedFile)
          {
              bedfile=champ.bedfile(resList_anno)
              #add pvalue to file name
              fileName1=paste(resultsDir,"/MVP_",adjPVal,"_",adjust.method,"adjust_",dim(resList)[1],".bed",sep="")
              write.table(bedfile,fileName1,row.names=F,col.names=F,quote=F, sep = "\t")
          }
          }
          resultsALL<- topTable(fit3, coef=1, number=dim(data)[1], adjust.method=adjust.method,p.value=1)
          resultsALL_anno <- data.frame(resultsALL,probe.features[match(resultsALL$ID, rownames(probe.features)),])
          
          data=data.frame(data)
          data<-cbind(data,C_AVG=apply(data[,which(colnames(data) %in% controls$Sample_Name)],1,mean))[1:nrow(data),]
          data<-cbind(data,T_AVG=apply(data[,which(colnames(data) %in% test$Sample_Name)],1,mean))[1:nrow(data),]
          data$deltaBeta<-abs(data$C_AVG-data$T_AVG)
          colnames(data)[length(data)-2]=paste(groupLabel[1],"_AVG",sep="")
          colnames(data)[length(data)-1]=paste(groupLabel[2],"_AVG",sep="")
          data1=data[(length(data)-2):length(data)]
          resultsALL_anno<-data.frame(resultsALL_anno,data1[match(rownames(data1),resultsALL_anno$ID),])
          colnames(resultsALL_anno)[1]="probeID"
          fileName2=paste(resultsDir,"/MVP_ALL_",groupLabel[1],"vs",groupLabel[2],"_",adjust.method,"adjust.txt",sep="")
          write.table(resultsALL_anno, fileName2 ,quote=F,sep="\t",row.names=F)
          return(results.file=resultsALL_anno)
      }
     else(print("There are no significant MVPs"))

      }
}

champ.bedfile<-function(info)
{
	bedfile=data.frame(matrix(NA,ncol=4,nrow=nrow(info)))

	#check for MAPINFO column
	colnames(bedfile)[1]="chr"
	colnames(bedfile)[2]="start"
	colnames(bedfile)[3]="end"
	colnames(bedfile)[4]="probe"
	bedfile$probe=rownames(info)
	bedfile$chr<-paste("chr",info$CHR,sep="")
	bedfile$start<-info$MAPINFO-1
	bedfile$end<-info$MAPINFO

	bedfile<-bedfile[order(bedfile$chr,bedfile$start),]
	return(bedfile)
}

###########DMR##############
champ.lasso <- function(fromFile=FALSE,uploadResults=FALSE,uploadFile="limma.txt",limma,beta.norm=norm$beta,pd=load$pd,filterXY=TRUE,image=TRUE,mafPol.lower=0, mafPol.upper=0.05,popPol="eur",lassoStyle="max",lassoRadius=2000, minSigProbesLasso=3, minDmrSep=1000, minDmrSize=0, adjPVal=0.05, adjust.method="BH",resultsDir=paste(getwd(),"resultsChamp",sep="/"),bedFile=FALSE,DMRpval=0.05,batchDone=FALSE,normSave)
{
	#probe.features<-NULL
	#probe.450K.VCs.af<-NULL

	message("Run DMR")
	if(uploadResults)
	{
		resultsFile <- read.table(uploadFile, header=T)
	}else if(fromFile){
		resultsFile<-limma
		
	}else
	{
		resultsFile=champ.MVP(beta.norm=beta.norm,pd=pd,adjPVal=adjPVal,bedFile=bedFile,adjust.method=adjust.method,resultsDir=resultsDir)
	}	
    myResults <- data.frame(resultsFile[,grep("adj.P.Val", colnames(resultsFile))])
    rownames(myResults) <- resultsFile$ID
    colnames(myResults) <- "adj.P.Val"
    myResults <- data.frame(myResults, probe.features[match(rownames(myResults), rownames(probe.features)),c(1:3,20)])

    if(dim(myResults)[1]==0)
    {
		message("Your dataset is empty and champ.lasso() cannot proceed")
		return(beta.norm)
	}else if(min(myResults$adj.P.Val)>adjPVal)
          {
            if(batchDone)
            {
                message("The adusted p-values in your dataset exceed the cutoff you have chosen champ.lasso() cannot proceed. champ.lasso will rerun with the normalised beta values without batchCorrect.")
                    champ.lasso(beta.norm=normSave,pd=pd,resultsDir=resultsDir,bedFile=bedFile,batchDone=F)
            }else{
                message("The adusted p-values in your dataset exceed the cutoff you have chosen and champ.lasso() cannot proceed.")
                return(beta.norm)
                }
        }else if(count(myResults$adj.P.Val<adjPVal)[which(count(myResults$adj.P.Val<adjPVal)$x==TRUE),][1,2] < 3){
		message("There are not enough MVPs in your dataset for champ.lasso() to proceed.")
		return(beta.norm)
	
            }else{
                #myResults <- data.frame(myResults, probe.features[match(rownames(myResults), rownames(probe.features)),c(1:3,20)])
    }
    
	if(filterXY) # converts $CHR into numeric for sorting purposes
	{
		myResults <- myResults[which(!(myResults$CHR=="X" | myResults$CHR=="Y")),]
		myResults$CHR <- as.numeric(as.character(myResults$CHR))
	}else
	{	
		myResults$CHR <- as.character(myResults$CHR)
		myResults$CHR <- replace(myResults$CHR, which(myResults$CHR == "X"), "23")
		myResults$CHR <- replace(myResults$CHR, which(myResults$CHR == "Y"), "24")
		myResults$CHR <- as.numeric(myResults$CHR)
	}

	pop.ref <- data.frame("popPol" = c("asn", "amr", "afr", "eur"), "col.f"=c(9, 13, 17, 21), "col.r"=c(10,14,18,22))
	vc <- data.frame("ID" = rownames(probe.450K.VCs.af), "pol.af.f" = probe.450K.VCs.af[,pop.ref$col.f[match(popPol, pop.ref$popPol)]], "pol.af.r" = probe.450K.VCs.af[,pop.ref$col.r[match(popPol, pop.ref$popPol)]])
	myPol <- cbind(vc$pol.af.f[match(rownames(myResults), vc$ID)], 1 - vc$pol.af.f[match(rownames(myResults), vc$ID)], vc$pol.af.r[match(rownames(myResults), vc$ID)], 1 - vc$pol.af.r[match(rownames(myResults), vc$ID)])
	myPol <- cbind(apply(myPol[, 1:2], 1, min), apply(myPol[, 3:4], 1, min))
	myResults <- data.frame(myResults, "pol.af.f" = myPol[,1], "pol.af.r" = myPol[,2])

	if(mafPol.lower == 0 & mafPol.upper == 0 )
	{
		myResults <- myResults[myResults$pol.af.f == 0 & myResults$pol.af.r == 0,]
	}else
	{
		myResults <- myResults[myResults$pol.af.f >= mafPol.lower &
		myResults$pol.af.f <= mafPol.upper |
		myResults$pol.af.r >= mafPol.lower &
		myResults$pol.af.r <= mafPol.upper,]
	}
	
	myResults <- myResults[order(myResults$CHR, myResults$MAPINFO),]
    
###Probe spacing and quantile derivation
	myResults.split <- split(myResults$MAPINFO, paste(myResults$CHR, myResults$arm), drop=T) #split by chromosome & arm
	myResults.split.rownames <- split(rownames(myResults), paste(myResults$CHR, myResults$arm), drop=T) # rownames split by chromosome
	probe.spacing <- lapply(myResults.split, diff)	# nearest probe calculations
	up.down <- data.frame("upstream" = unlist(lapply(probe.spacing, function(x) c(0, x))),
    "downstream" = unlist(lapply(probe.spacing, function(x) c(x,0))), row.names = unlist(myResults.split.rownames))
	up.down$nrst.probe <- apply(up.down, 1, min)
	rpl <- apply(up.down[which(up.down$nrst.probe == 0),],1,max)
	up.down$nrst.probe <- replace(up.down$nrst.probe, which(up.down$nrst.probe == 0), rpl)
	myResults <- data.frame(myResults, "nrst.probe" = up.down$nrst.probe[match(rownames(myResults), rownames(up.down))])

    
	lasso.quantiles <- do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.rel), function(x) ecdf(x)(lassoRadius)))
	if(lassoStyle == "max")
	{
		value.lasso.quantile <- min(lasso.quantiles)
	}else
	{
		value.lasso.quantile <- max(lasso.quantiles)
	}
    
	lasso.radii <- round(do.call(rbind, lapply(split(myResults$nrst.probe, myResults$feat.rel), function(x) quantile(x, value.lasso.quantile, na.rm=T))))
	myResults$lasso.radii <- lasso.radii[match(myResults$feat.rel, rownames(lasso.radii))]
	
	myResults.sig <- myResults[myResults$adj.P.Val < adjPVal,]
	myResults.sig.splitarm <- split(myResults.sig[,c(3,ncol(myResults.sig))], paste(myResults.sig$CHR, myResults.sig$arm), drop=T) #split MAPINFO & lasso.radii by chromosome & arm
	myResults.sig.rownames <- split(rownames(myResults.sig), paste(myResults.sig$CHR, myResults.sig$arm), drop=T) # rownames split by chromosome

	no.pr <- vector("list", length(myResults.sig.splitarm)) # counts no. of probes captured by lasso
	for (k in 1:length(myResults.sig.splitarm))
	{
		for (i in 1:nrow(myResults.sig.splitarm[[k]]))
		{
			no.pr[[k]][i] <- 	length(which(myResults.sig.splitarm[[k]][,1] > myResults.sig.splitarm[[k]][i,1] - myResults.sig.splitarm[[k]][i,2] &
            myResults.sig.splitarm[[k]][,1] < myResults.sig.splitarm[[k]][i,1] + myResults.sig.splitarm[[k]][i,2]))
		}
	}
    
        
	tmp1 <- data.frame("lassoRadius" = do.call("rbind", myResults.sig.splitarm)[,2], "lasso.captured.probes"=unlist(no.pr), row.names = unlist(myResults.sig.rownames))
    
    if(minSigProbesLasso > max(tmp1$lasso.captured.probes))
    {
        message("Your lassos have not managed to capture at least ", minSigProbesLasso, " probes. No DMRs can be found")
        return(beta.norm)
    } else{
    
    tmp2 <- data.frame(myResults.sig[match(rownames(tmp1), rownames(myResults.sig)), ], "lasso.captured.probes"=tmp1$lasso.captured.probes)
    tmp3 <- tmp2[tmp2$lasso.captured.probes >= minSigProbesLasso,] # retains rows (probes) that capture no less than "no.probes"
    tmp3 <- tmp3[order(tmp3$CHR, tmp3$MAPINFO),] # orders file by chromosome then position
	
    lasso.coord <- data.frame("CHR" = tmp3$CHR, "lasso.start"=tmp3$MAPINFO - tmp3$lasso.radii, "lasso.end"=tmp3$MAPINFO + tmp3$lasso.radii)
    lasso.coord <- lasso.coord[order(lasso.coord$CHR, lasso.coord$lasso.start),]
	
    lasso.diffs <- numeric(nrow(lasso.coord)) # looking for overlapping lassos
    for (i in 2:length(lasso.diffs))
    {
        lasso.diffs[i] <- lasso.coord$lasso.start[i] - lasso.coord$lasso.end[i-1]
    }
        
    dmr.marker <- numeric(length(lasso.diffs))
    for (i in 2:length(dmr.marker))
    {
        dmr.marker[i] <- ifelse(tmp3$CHR[i] != tmp3$CHR[i-1] | lasso.diffs[i] > minDmrSep, 0, 1)
    }
	
    dmr.no <- numeric(length(dmr.marker))
    dmr.no[1] <- 1
    for (i in 2:length(dmr.marker))
    {
        dmr.no[i] <- ifelse(dmr.marker[i] == 0, dmr.no[i-1]+1, dmr.no[i-1])
    }
	
    dmrBylasso <- data.frame(lasso.coord, dmr.no)
    dmrs <- data.frame("chr"=unlist(lapply(split(dmrBylasso$CHR,dmrBylasso$dmr.no),unique)), "dmr.no"=unique(dmrBylasso$dmr.no), "dmr.start"=unlist(lapply(split(dmrBylasso[,2:3], dmrBylasso$dmr.no),min)), "dmr.end"=unlist(lapply(split(dmrBylasso[,2:3], dmrBylasso$dmr.no),max)))
    dmrs$dmr.size <- dmrs$dmr.end - dmrs$dmr.start +1
    dmrs <- dmrs[dmrs$dmr.size >= minDmrSize,]
    rownames(dmrs) <- 1:nrow(dmrs) # renames rows after removing small DMRs
    dmrs$dmr.no <- 1:nrow(dmrs) # renames DMRs after removing small DMRs
        
        dmr.probes <- vector("list", nrow(dmrs)) # collects all probes (sig and nonsig) within DMR boundaries
        names(dmr.probes) <- dmrs$dmr.no
        for (i in 1:length(dmr.probes))
        {
            dmr.probes[[i]] <- which(myResults$CHR == dmrs$chr[i] & myResults$MAPINFO >= dmrs$dmr.start[i] & myResults$MAPINFO <= dmrs$dmr.end[i])
        }
        dmr.probes.lengths <- unlist(lapply(dmr.probes, length))
        dmr.probes <- data.frame("dmr.no"=rep(as.numeric(names(dmr.probes)), dmr.probes.lengths), "probe.no"=unlist(dmr.probes))
        
        tmp4 <- data.frame("ID" = rownames(myResults)[dmr.probes$probe.no], myResults[dmr.probes$probe.no,], dmrs[dmr.probes$dmr.no,-1], row.names=1:nrow(dmr.probes))
        
        if(filterXY) # retro-converts $CHR into a factor
        {
            tmp4$CHR <- as.factor(as.character(tmp4$CHR))
        }else
        {
            tmp4$CHR <- as.character(tmp4$CHR)
            tmp4$CHR <- replace(tmp4$CHR, which(tmp4$CHR == 23), "X")
            tmp4$CHR <- replace(tmp4$CHR, which(tmp4$CHR == 24), "Y")
            tmp4$CHR <- as.factor(tmp4$CHR)
        }
	
        dmr.beta <- split(as.data.frame(beta.norm[match(tmp4$ID, rownames(beta.norm)),]), tmp4$dmr.no)
		corel <- lapply(dmr.beta, function(x) cor(t(x)))
		weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))
		dmr.ind.p <- split(tmp4$adj.P.Val, tmp4$dmr.no)
		dmr.qp <- lapply(dmr.ind.p, qnorm)
		dmr.qp.w <- mapply("*", dmr.qp, weights)
		dmr.stat <- lapply(dmr.qp.w, sum)
		dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
		dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
		dmr.rep <- as.numeric(summary(dmr.ind.p)[,1])
		dmr.p <- dmr.p*length(dmr.p)/rank(dmr.p)
        dmr.p <- rep(dmr.p, dmr.rep)
    
    if(image)
	{
		plotThreshold(myResults,lasso.radii,value.lasso.quantile,lassoRadius,lassoStyle,resultsDir)
	}
		
    message("You have found ",max(dmr.probes$dmr.no)," significant DMRs.")
	dmrList=data.frame(tmp4, "dmr.p" = dmr.p)
	colnames(dmrList)[1]="probeID"
	fileName=paste(resultsDir,"/DMR_",adjPVal,".txt",sep="")
	write.table(dmrList, fileName ,quote=F,sep="\t",row.names=F)
	
	########make bedfile
	#lasso_test=data.frame(lasso$CHR,lasso$dmr.start,lasso$dmr.stop,lasso$dmr.no)
	#lasso_test2=unique(lasso_test)
	#colnames(lasso_test2)<-c("chr","start","end","dmr")
	#lasso_test2$chr<-paste("chr",lasso_test2$chr,sep="")
	#write.table(lasso_test2,"pr_DMR.bed",row.names=F,col.names=F,quote=F, sep = "\t")
	
	return(dmrList)
    }
}

########################################### PLOT TO DISCOVER THRESHOLD RADII ### (see #note above)
plotThreshold <- function(myResults,lasso.radii,value.lasso.quantile,lassoRadius,lassoStyle,resultsDir)
{
#probe.features<-NULL
#myResults<-NULL

	no.feature <- length(summary(as.factor(probe.features$feature.1)))
	no.relation <- length(summary(as.factor(probe.features$RELATION_TO_UCSC_CPG_ISLAND))) - 2 # merges North- & South- shore & shelf relations
	no.feat.rel <- length(summary(as.factor(probe.features$feat.rel)))
    
	spc.thrd.qt <- matrix(nrow=99, ncol=no.feat.rel)
	for (i in 1:99)
	{
		spc.thrd.qt[i,] <- as.numeric(tapply(myResults$nrst.probe, myResults$feat.rel, function(x) quantile(x,i/100, na.rm=T))) # calculates nearest probe split by feature
	}
    
    imageName <- paste(resultsDir,"myLassos.pdf",sep="/")
    pdf(imageName,width=9,height=9)
	par(bty="n", adj=0.5, mar=c(8,4,7,2)+0.1)
	boxplot(myResults$nrst.probe~myResults$feat.rel, outline=F, ylab="nearest neighbour [bp]", col=rep(rainbow(no.feature),each=no.relation), las=2, lwd=0.2)
    axis(side=3, at=c(1:no.feat.rel), labels=paste("N =", format(as.numeric(tapply(myResults$nrst, myResults$feat.rel, length)), big.mark=",", scientific=F), sep=""), tick=FALSE, cex.axis=0.8, las=2)
    title(main="a", adj=0, cex.main=2)
	plot(log10(spc.thrd.qt[,1]), ylim=log10(range(spc.thrd.qt)), xlab="quantile", ylab="lasso radius [bp]", lty=1, lwd=2, col=rainbow(no.feature)[1], type="l", xaxt="n", yaxt="n")
    axis(side=1, at=seq(0,100,10), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,5,1), labels=c(0,10,100,1000,10000,100000), las=2)
    for (i in 2:no.feat.rel)
    {
        lines(log10(spc.thrd.qt[,i]), lty=rep(c(1,2,3,6),no.feature)[i], lwd=2, col=rep(rainbow(no.feature),each=no.relation)[i])
    }
    legend(0,5,legend=c(unique(sapply(strsplit(names(summary(as.factor(probe.features$feat.rel))), "_"), "[[",1)),unique(sapply(strsplit(names(summary(as.factor(probe.features$feat.rel))), "_"), "[[",2))), col=c(rainbow(7),rep("black",4)), lty=c(rep(1,7),c(1,2,3,6)), lwd=2, bty="n", cex=0.8)
    if(lassoStyle == "max")
    {
        segments(c(0,value.lasso.quantile*100), rep(log10(lassoRadius), 2), rep(value.lasso.quantile*100,2), c(log10(lassoRadius),0), lty=2)
    }else
    {
        segments(c(0,value.lasso.quantile*100), c(log10(lassoRadius), log10(max(lasso.radii))), rep(value.lasso.quantile*100,2), c(log10(lassoRadius),0), lty=2)
    }
    title(main="b", adj=0, cex.main=2)
	plot(c(1,28), y=c(range(0.3*sqrt(lasso.radii))[1]*0.8, range(0.3*sqrt(lasso.radii))[2]*1.2), type="n", xaxt="n", xlab="", yaxt="n", ylab="lasso radius [bp]", main=paste("lasso quantile = ", round(value.lasso.quantile,2), sep=""), bty="n")
    segments(1:28, rep(0,28), 1:28, 0.3*sqrt(lasso.radii), lty=3, col="grey")
    points(1:28, 0.3*sqrt(lasso.radii), pch=16, cex=0.3*sqrt(lasso.radii), col=rep(rainbow(7,alpha=0.5), each=4))
    text(1:28, 0.3*sqrt(lasso.radii), lasso.radii, pos=3, cex=0.8)
    axis(1, at=1:28,rownames(lasso.radii), las=2, cex.axis=0.8, tick=F)
    axis(2, at=c(0,max(0.3*sqrt(lasso.radii))), labels=F)
    title(main="c", adj=0, cex.main=2)
	dev.off()
}


###########CNA##############
champ.CNA<-function(intensity=load$intensity, pd=load$pd, loadFile=F, batchCorrect=T, file="intensity.txt", control=FALSE,resultsDir=paste(getwd(),"resultsChamp",sep="/"),plotCNA=T) 
{
	#probe.features<-NULL
	#normalize.quantiles<-NULL
    #rm(normalize.quantiles)
	#control.intsqnlog<-NULL
	CNA<-NA
	rm(CNA)
	smooth.CNA<-NA
	rm(smooth.CNA)
	segment<-NA
	rm(segment)
	
	if(require(preprocessCore)){if(require(DNAcopy)){
	message("Run CNA")
	newDir=paste(resultsDir,"CNA",sep="/")
	if(plotCNA){if(!file.exists(resultsDir)){dir.create(resultsDir)}
	if(!file.exists(newDir))
	{
		dir.create(newDir)
	}}
		
	if(loadFile)
	{
		ints = read.table(file,row.names =T, sep = "\t")
	}else
	{
		ints=intensity
	}
	
	#Extracts names of samples 
	names<-colnames(ints)

	#Quantile normalises intensities	
	intsqn<-normalize.quantiles(as.matrix(ints))
	colnames(intsqn)<-names
	
	if(batchCorrect)
	{
		message("Run batch correction for CNA")
		combat=champ.runCombat(beta.c=intsqn,pd=pd,logitTrans=FALSE)
		intsqn=combat$beta
	}

	#Calculates Log2 
	intsqnlog<-log2(intsqn)
	
	if(control)
	{
		#separates case from control(reference sample/samples)
		#case.intsqnlog<-intsqnlog[,1:length(names)-1]
		#control.intsqnlog<-intsqnlog[,length(names)]
		message("This option is not yet available")
		control=F
			
	}else
	{
		#Creates alternate reference sample from rowMeans if proper reference /control is not available 
		case.intsqnlog<-intsqnlog[,1:length(names)]
		ref.intsqnlog<-rowMeans(intsqnlog)
	}
	#Generates Log2Ratio for case v control/reference
	intsqnlogratio<-intsqnlog
	colnames(intsqnlogratio)<-names
	for(i in 1:ncol(case.intsqnlog))
	{
		if(control)
		{
			intsqnlogratio[,i]<-case.intsqnlog[,i]-control.intsqnlog
		}else 
		{
			intsqnlogratio[,i]<-case.intsqnlog[,i]-ref.intsqnlog
		}
	}

	ints <- data.frame(ints, probe.features$MAPINFO[match(rownames(ints), rownames(probe.features))])
	names(ints)[length(ints)] <- "MAPINFO"
	ints <- data.frame(ints, probe.features$CHR[match(rownames(ints), rownames(probe.features))])
	names(ints)[length(ints)] <- "CHR"

	#Replaces Chr X and Y with 23 and 24
	levels(ints$CHR)[levels(ints$CHR)=='X']='23'
	levels(ints$CHR)[levels(ints$CHR)=='Y']='24'

	#converts CHR factors to numeric 
	CHR<-as.numeric(levels(ints$CHR))[ints$CHR]

	#need to copy in MAPINFO
	MAPINFO<-ints$MAPINFO
	
	#Runs CNA and generates individual DNA Copy Number profiles
	for(i in 1:ncol(intsqnlogratio))
	{
		CNA.object <- CNA(cbind(intsqnlogratio[,i]), CHR, MAPINFO ,data.type = "logratio", sampleid = paste(names[i],"_qn"))
		smoothed.CNA.object <- smooth.CNA(CNA.object)
		segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2)
		if(plotCNA)
		{
			imageName<-paste(names[i],"_qn.jpg",sep="")
			imageName=paste(newDir,imageName,sep="/")
			jpeg(imageName)
			plot(segment.smoothed.CNA.object, plot.type = "w", ylim=c(-6,6))
			dev.off()
		}
		seg<-segment.smoothed.CNA.object$output
        
		table_name<-paste(newDir,"/",names[i],"_qn.txt",sep="")
		write.table(seg,table_name, sep="\t", col.names=T, row.names=F, quote=FALSE)
	}}}
}

##############################################################
#######################BEGIN SCRIPTS##########################
##############################################################

### BMIQ.R & CheckBMIQ.R
### This function adjusts for the type-2 bias in Illumina Infinium 450k data.
### Author: Andrew Teschendorff
### Date v_1.1: Nov 2012
### Date v_1.2: 6th Apr 2013
### Date v_1.3: 29th May 2013

### SUMMARY
### BMIQ is an intra-sample normalization procedure, adjusting for the bias in type-2 probe values, using a 3-step procedure published in Teschendorff AE et al "A Beta-Mixture Quantile Normalization method for correcting probe design bias in Illumina Infinium 450k DNA methylation data", Bioinformatics 2012 Nov 21.


### INPUT:
### beta.v: vector consisting of beta-values for a given sample. NAs are not allowed, so these must be removed or imputed prior to running BMIQ. Beta-values that are exactly 0 or 1 will be replaced by the minimum positive above zero or maximum value below 1, respectively.
### design.v: corresponding vector specifying probe design type (1=type1,2=type2). This must be of the same length as beta.v and in the same order.
### doH: perform normalization for hemimethylated type2 probes. By default TRUE.
### nfit: number of probes of a given design to use for the fitting. Default is 10000. Smaller values will make BMIQ run faster at the expense of a small loss in accuracy. For most applications, even 5000 is ok.
### nL: number of states in beta mixture model. 3 by default. At present BMIQ only works for nL=3.
### th1.v: thresholds used for the initialisation of the EM-algorithm, they should represent buest guesses for calling type1 probes hemi-methylated and methylated, and will be refined by the EM algorithm. Default values work well in most cases.
### th2.v: thresholds used for the initialisation of the EM-algorithm, they should represent buest guesses for calling type2 probes hemi-methylated and methylated, and will be refined by the EM algorithm. By default this is null, and the thresholds are estimated based on th1.v and a modified PBC correction method.
### niter: maximum number of EM iterations to do. By default 5.
### tol: tolerance threshold for EM algorithm. By default 0.001.
### plots: logical specifying whether to plot the fits and normalised profiles out. By default TRUE.
### sampleID: the ID of the sample being normalised.

### OUTPUT
### A list with the following elements:
### nbeta: the normalised beta-profile for the sample
### class1: the assigned methylation state of type1 probes
### class2: the assigned methylation state of type2 probes
### av1: mean beta-values for the nL classes for type1 probes.
### av2: mean beta-values for the nL classes for type2 probes.
### hf: the "Hubble" dilation factor
### th1: estimated thresholds used for type1 probes
### th2: estimated thresholds used for type2 probes

require(RPMM);

BMIQ <- function(beta.v,design.v,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1){

type1.idx <- which(design.v==1);
type2.idx <- which(design.v==2);

beta1.v <- beta.v[type1.idx];
beta2.v <- beta.v[type2.idx];

### check if there are exact 0's or 1's. If so, regularise using minimum positive and maximum below 1 values.
if(min(beta1.v)==0){
  beta1.v[beta1.v==0] <- min(setdiff(beta1.v,0));
}
if(min(beta2.v)==0){
  beta2.v[beta2.v==0] <- min(setdiff(beta2.v,0));
}
if(max(beta1.v)==1){
  beta1.v[beta1.v==1] <- max(setdiff(beta1.v,1));
}
if(max(beta2.v)==1){
  beta2.v[beta2.v==1] <- max(setdiff(beta2.v,1));
}

### estimate initial weight matrix from type1 distribution
w0.m <- matrix(0,nrow=length(beta1.v),ncol=nL);
w0.m[which(beta1.v <= th1.v[1]),1] <- 1;
w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] <- 1;
w0.m[which(beta1.v > th1.v[2]),3] <- 1;

### fit type1
print("Fitting EM beta mixture to type1 probes");
rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
em1.o <- blc(matrix(beta1.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
subsetclass1.v <- apply(em1.o$w,1,which.max);
subsetth1.v <- c(mean(c(max(beta1.v[rand.idx[subsetclass1.v==1]]),min(beta1.v[rand.idx[subsetclass1.v==2]]))),mean(c(max(beta1.v[rand.idx[subsetclass1.v==2]]),min(beta1.v[rand.idx[subsetclass1.v==3]]))));
class1.v <- rep(2,length(beta1.v));
class1.v[which(beta1.v < subsetth1.v[1])] <- 1;
class1.v[which(beta1.v > subsetth1.v[2])] <- 3;
nth1.v <- subsetth1.v;
print("Done");

### generate plot from estimated mixture
if(plots){
print("Check");
tmpL.v <- as.vector(rmultinom(1:nL,length(beta1.v),prob=em1.o$eta));
tmpB.v <- vector();
for(l in 1:nL){
  tmpB.v <- c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
}

pdf(paste("Type1fit-",sampleID,".pdf",sep=""),width=6,height=4);
plot(density(beta1.v));
d.o <- density(tmpB.v);
points(d.o$x,d.o$y,col="green",type="l")
legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
dev.off();
}



### Estimate Modes 
d1U.o <- density(beta1.v[class1.v==1])
d1M.o <- density(beta1.v[class1.v==3])
mod1U <- d1U.o$x[which.max(d1U.o$y)]
mod1M <- d1M.o$x[which.max(d1M.o$y)]
d2U.o <- density(beta2.v[which(beta2.v<0.4)]);
d2M.o <- density(beta2.v[which(beta2.v>0.6)]);
mod2U <- d2U.o$x[which.max(d2U.o$y)]
mod2M <- d2M.o$x[which.max(d2M.o$y)]


### now deal with type2 fit
th2.v <- vector();
th2.v[1] <- nth1.v[1] + (mod2U-mod1U);
th2.v[2] <- nth1.v[2] + (mod2M-mod1M);

### estimate initial weight matrix 
w0.m <- matrix(0,nrow=length(beta2.v),ncol=nL);
w0.m[which(beta2.v <= th2.v[1]),1] <- 1;
w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] <- 1;
w0.m[which(beta2.v > th2.v[2]),3] <- 1;

print("Fitting EM beta mixture to type2 probes");
rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
em2.o <- blc(matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
print("Done");

### for type II probes assign to state (unmethylated, hemi or full methylation)
subsetclass2.v <- apply(em2.o$w,1,which.max);
subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v==1]]),min(beta2.v[rand.idx[subsetclass2.v==2]])),mean(max(beta2.v[rand.idx[subsetclass2.v==2]]),min(beta2.v[rand.idx[subsetclass2.v==3]])));
class2.v <- rep(2,length(beta2.v));
class2.v[which(beta2.v < subsetth2.v[1])] <- 1;
class2.v[which(beta2.v > subsetth2.v[2])] <- 3;


### generate plot
if(plots){
tmpL.v <- as.vector(rmultinom(1:nL,length(beta2.v),prob=em2.o$eta));
tmpB.v <- vector();
for(lt in 1:nL){
  tmpB.v <- c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
}
pdf(paste("Type2fit-",sampleID,".pdf",sep=""),width=6,height=4);
plot(density(beta2.v));
d.o <- density(tmpB.v);
points(d.o$x,d.o$y,col="green",type="l")
legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
dev.off();
}

classAV1.v <- vector();classAV2.v <- vector();
for(l in 1:nL){
  classAV1.v[l] <-  em1.o$mu[l,1];
  classAV2.v[l] <-  em2.o$mu[l,1];
}

### start normalising type2 probes
print("Start normalising type 2 probes");
nbeta2.v <- beta2.v;
### select U probes
lt <- 1;
selU.idx <- which(class2.v==lt);
selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
### find prob according to typeII distribution
p.v <- pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
### find corresponding quantile in type I distribution
q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
nbeta2.v[selUR.idx] <- q.v;
p.v <- pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
### find corresponding quantile in type I distribution
q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
nbeta2.v[selUL.idx] <- q.v;

### select M probes
lt <- 3;
selM.idx <- which(class2.v==lt);
selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
### find prob according to typeII distribution
p.v <- pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
### find corresponding quantile in type I distribution
q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
nbeta2.v[selMR.idx] <- q.v;


if(doH){ ### if TRUE also correct type2 hemimethylated probes
### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
lt <- 2;
selH.idx <- c(which(class2.v==lt),selML.idx);
minH <- min(beta2.v[selH.idx])
maxH <- max(beta2.v[selH.idx])
deltaH <- maxH - minH;
#### need to do some patching
deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])

## new maximum of H probes should be
nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM;
## new minimum of H probes should be
nminH <- max(nbeta2.v[selU.idx]) + deltaUH;
ndeltaH <- nmaxH - nminH;

### perform conformal transformation (shift+dilation)
## new_beta_H(i) = a + hf*(beta_H(i)-minH);
hf <- ndeltaH/deltaH ;
### fix lower point first
nbeta2.v[selH.idx] <- nminH + hf*(beta2.v[selH.idx]-minH);

}

pnbeta.v <- beta.v;
pnbeta.v[type1.idx] <- beta1.v;
pnbeta.v[type2.idx] <- nbeta2.v;

### generate final plot to check normalization
if(plots){
 print("Generating final plot");
 d1.o <- density(beta1.v);
 d2.o <- density(beta2.v);
 d2n.o <- density(nbeta2.v);
 ymax <- max(d2.o$y,d1.o$y,d2n.o$y);
 pdf(paste("CheckBMIQ-",sampleID,".pdf",sep=""),width=6,height=4)
 plot(density(beta2.v),type="l",ylim=c(0,ymax),xlim=c(0,1));
 points(d1.o$x,d1.o$y,col="red",type="l");
 points(d2n.o$x,d2n.o$y,col="blue",type="l");
 legend(x=0.5,y=ymax,legend=c("type1","type2","type2-BMIQ"),bty="n",fill=c("red","black","blue"));
 dev.off();
}

print(paste("Finished for sample ",sampleID,sep=""));

return(list(nbeta=pnbeta.v,class1=class1.v,class2=class2.v,av1=classAV1.v,av2=classAV2.v,hf=hf,th1=nth1.v,th2=th2.v));

}



CheckBMIQ <- function(beta.v,design.v,pnbeta.v){### pnbeta is BMIQ normalised profile

type1.idx <- which(design.v==1);
type2.idx <- which(design.v==2);

beta1.v <- beta.v[type1.idx];
beta2.v <- beta.v[type2.idx];
pnbeta2.v <- pnbeta.v[type2.idx];
  

}

###########PBC#################
DoPBC <- function(beta.m,design.v){

  mval.m <- log2(beta.m/(1-beta.m));
  type1.idx <- which(design.v==1)
  type2.idx <- which(design.v==2)  
  mvalT.m <- mval.m;
  for(s in 1:ncol(beta.m)){

    neg.idx <- which(mval.m[,s]<0);
    pos.idx <- which(mval.m[,s]>0);

    neg1.idx <- intersect(neg.idx,type1.idx);
    neg2.idx <- intersect(neg.idx,type2.idx);

    pos1.idx <- intersect(pos.idx,type1.idx);
    pos2.idx <- intersect(pos.idx,type2.idx);

    d.o <- density(mval.m[neg2.idx,s],kernel="gaussian",bw=0.5);
    peakU2 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[pos2.idx,s],kernel="gaussian",bw=0.5);
    peakM2 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[neg1.idx,s],kernel="gaussian",bw=0.5);
    peakU1 <- abs(d.o$x[which.max(d.o$y)]);

    d.o <- density(mval.m[pos1.idx,s],kernel="gaussian",bw=0.5);
    peakM1 <- abs(d.o$x[which.max(d.o$y)]);

    mvalT.m[neg2.idx,s] <- (mval.m[neg2.idx,s]/peakU2)*peakU1;
    mvalT.m[pos2.idx,s] <- (mval.m[pos2.idx,s]/peakM2)*peakM1;
    print(paste("Done for sample ",s,sep=""));
  }

  betaT.m <- 2^mvalT.m/(2^mvalT.m+1);
  return(betaT.m);
}

### EstDimRMTv2.R

thdens <- function(Q,sigma2,ns){

  lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
  lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
  
  delta <- lambdaMAX - lambdaMIN;#  print(delta);
  
  roundN <- 3;
  step <- round(delta/ns,roundN);
  while(step==0){
    roundN <- roundN+1;
    step <- round(delta/ns,roundN);
  }
  
  lambda.v <- seq(lambdaMIN,lambdaMAX,by=step);
  dens.v <- vector();
  ii <- 1;
  for(i in lambda.v){
    dens.v[ii] <- (Q/(2*pi*sigma2))*sqrt( (lambdaMAX-i)*(i-lambdaMIN) )/i;
    ii <- ii+1;
  }

  return(list(min=lambdaMIN,max=lambdaMAX,step=step,lambda=lambda.v,dens=dens.v));
}

GenPlot <- function(thdens.o,estdens.o,evalues.v){


minx <- min(min(thdens.o$lambda),min(evalues.v));
maxx <- max(max(thdens.o$lambda),max(evalues.v));
miny <- min(min(thdens.o$dens),min(estdens.o$y));
maxy <- max(max(thdens.o$dens),max(estdens.o$y));


#plot(thdens.o$lambda,thdens.o$dens,xlim=c(minx,maxx),ylim=c(miny,maxy),type="b",col="green");
#i <- min(which(estdens.o$x > min(evalues.v)));
#f <- max(which(estdens.o$x < max(evalues.v)));
#points(x=estdens.o$x[i:f],y=estdens.o$y[i:f],type="b",col="red");

}


EstDimRMTv2 <- function(data.m){

 ### standardise matrix
 M <- data.m;
 for(c in 1:ncol(M)){
  M[,c] <- (data.m[,c]-mean(data.m[,c]))/sqrt(var(data.m[,c]));
 }
 sigma2 <- var(as.vector(M));
 Q <- nrow(data.m)/ncol(data.m);
 thdens.o <- thdens(Q,sigma2,ncol(data.m));

 C <- 1/nrow(M) * t(M) %*% M;

 eigen.o <- eigen(C,symmetric=TRUE);
 estdens.o <- density(eigen.o$values,from=min(eigen.o$values),to=max(eigen.o$values),cut=0);

 GenPlot(thdens.o,estdens.o,eigen.o$values);
 intdim <- length(which(eigen.o$values > thdens.o$max));

 return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o));
 
}

champ.ComBat <- function (dat, batch, mod, numCovs = NULL, par.prior = TRUE,
prior.plots = FALSE)
{
    mod = cbind(mod, batch)
    check = apply(mod, 2, function(x) all(x == 1))
    mod = as.matrix(mod[, !check])
    colnames(mod)[ncol(mod)] = "Batch"
    if (sum(check) > 0 & !is.null(numCovs))
    numCovs = numCovs - 1
    design <- design.mat(mod, numCov = numCovs)
    batches <- list.batch(mod)
    n.batch <- length(batches)
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    NAs = any(is.na(dat))
    if (NAs) {
        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"),
        sep = " ")
    }
    cat("Standardizing Data across genes\n")

    tryCatch(
    
    if (!NAs) {                
        B.hat <- solve(t(design) %*% design) %*% t(design) %*%
        t(as.matrix(dat))
    }
    else {
        B.hat = apply(dat, 1, Beta.NA, design)
    },
    error=function(e)
    {
        message("Combat failed...Your slides may be confounded with batch or with each other. Analysis will proceed without batch correction\n")
        
    })
    if(exists("B.hat"))
    {
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array,
        n.array)
    }
    else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var,
        na.rm = T)
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1,
    n.array)))
    cat("Fitting L/S model and finding priors\n")
    batch.design <- design[, 1:n.batch]
    if (!NAs) {
        gamma.hat <- solve(t(batch.design) %*% batch.design) %*%
        t(batch.design) %*% t(as.matrix(s.data))
    }
    else {
        gamma.hat = apply(s.data, 1, Beta.NA, batch.design)
    }
    delta.hat <- NULL
    for (i in batches) {
        delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var,
        na.rm = T))
    }
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, aprior)
    b.prior <- apply(delta.hat, 1, bprior)
    if (prior.plots & par.prior) {
        par(mfrow = c(2, 2))
        tmp <- density(gamma.hat[1, ])
        plot(tmp, type = "l", main = "Density Plot")
        xx <- seq(min(tmp$x), max(tmp$x), length = 100)
        lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
        qqnorm(gamma.hat[1, ])
        qqline(gamma.hat[1, ], col = 2)
        tmp <- density(delta.hat[1, ])
        invgam <- 1/rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
        tmp1 <- density(invgam)
        plot(tmp, typ = "l", main = "Density Plot", ylim = c(0,
        max(tmp$y, tmp1$y)))
        lines(tmp1, col = 2)
        qqplot(delta.hat[1, ], invgam, xlab = "Sample Quantiles",
        ylab = "Theoretical Quantiles")
        lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
        title("Q-Q Plot")
    }
    gamma.star <- delta.star <- NULL
    if (par.prior) {
        cat("Finding parametric adjustments\n")
        for (i in 1:n.batch) {
            temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i,
            ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
            b.prior[i])
            gamma.star <- rbind(gamma.star, temp[1, ])
            delta.star <- rbind(delta.star, temp[2, ])
        }
    }
    else {
        cat("Finding nonparametric adjustments\n")
        for (i in 1:n.batch) {
            temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
            gamma.hat[i, ], delta.hat[i, ])
            gamma.star <- rbind(gamma.star, temp[1, ])
            delta.star <- rbind(delta.star, temp[2, ])
        }
    }
    cat("Adjusting the Data\n")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,
        ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1,
        n.batches[j])))
        j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1,
    n.array)))) + stand.mean
    return(bayesdata)
    }
    else(return())

}

design.mat <- function(mod, numCov){
    
	tmp <- which(colnames(mod) == 'Batch')
	tmp1 <- as.factor(mod[,tmp])
	cat("Found",nlevels(tmp1),'batches\n')
	design <- build.design(tmp1,start=1)
    
	if(!is.null(numCov)) {
		theNumCov = as.matrix(mod[,numCov])
		mod0 = as.matrix(mod[,-c(numCov,tmp)])
	} else 	mod0 = as.matrix(mod[,-tmp])
    
	ncov <- ncol(mod0)
    
	cat("Found",ncov,' categorical covariate(s)\n')
	if(!is.null(numCov)) cat("Found",ncol(theNumCov),' continuous covariate(s)\n')
	if(ncov>0){
		for (j in 1:ncov){
			tmp1 <- as.factor(as.matrix(mod0)[,j])
			design <- build.design(tmp1,des=design)
        }
    }
	if(!is.null(numCov)) design = cbind(design,theNumCov)
	return(design)
    
}
build.design <- function(vec, des=NULL, start=2){
	tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
	for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
	cbind(des,tmp)
}

list.batch <- function(mod){
	tmp1 <- as.factor(mod[,which(colnames(mod) == 'Batch')])
	batches <- NULL
	for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
	batches
}

aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}

bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}
it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
	n <- apply(!is.na(sdat),1,sum)
	g.old <- g.hat
	d.old <- d.hat
	change <- 1
	count <- 0
	while(change>conv){
		g.new <- postmean(g.hat,g.bar,n,d.old,t2)
		sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
		d.new <- postvar(sum2,n,a,b)
		change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
		g.old <- g.new
		d.old <- d.new
		count <- count+1
    }
	#cat("This batch took", count, "iterations until convergence\n")
	adjust <- rbind(g.new, d.new)
	rownames(adjust) <- c("g.star","d.star")
	adjust
}
