# use latest champ in lastest R version. 

DMR.CHAMP<-function(){
  library("ChAMP")
  Rlt<-list()
  setwd("/home/shg047/bak/zh/data/zhu2016/dataset/idat")
  samplesheetFile=list.files(pattern="*csv")
  sampleID=unlist(strsplit(samplesheetFile,"[.]"))[1]
  # be sure to remove header of the samplesheet
  sampleinfo<-read.table(samplesheetFile,head=T,sep="\t",as.is=T)
  sampleinfo
  system("mkdir resultsChamp")
  # champ.process(fromIDAT = TRUE, fromFile = FALSE, directory = getwd(), resultsDir =paste(getwd(), "resultsChamp", sep = "/"), methValue = "B", filterDetP = TRUE,
  #              detPcut = 0.01, filterXY = TRUE, removeDetP = 0, filterBeads = TRUE, beadCutoff =
  #                0.05, filterNoCG = FALSE, QCimages = TRUE, batchCorrect = FALSE, runSVD =
  #                TRUE, studyInfo = FALSE, infoFactor = c(), norm = "BMIQ", adjust.method = "BH",
  #              adjPVal = 0.05, runDMR = TRUE, runCNA = FALSE, plotBMIQ = FALSE, DMRpval = 0.05,
  #              sampleCNA=FALSE,plotSample = TRUE,groupFreqPlots=TRUE,freqThreshold=0.3, bedFile
  #              = TRUE, methProfile = TRUE, controlProfile = FALSE)
  myLoad<-champ.load(directory = getwd(),filterBeads=TRUE,QCimages = F)
  # champ.CNA(intensity = myLoad$intensity, pd = myLoad$pd, loadFile = FALSE, batchCorrect = TRUE,
  #           file = "intensity.txt", resultsDir = paste(getwd(), "resultsChamp", sep = "/"),
  #           sampleCNA=TRUE, plotSample=TRUE, filterXY = TRUE, groupFreqPlots=TRUE,freqThreshold=0.3,
  #           control=TRUE,controlGroup="Control")
  # turn off the plot if you run it on the server since of the problem of X11
  myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = samplesheet, resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = TRUE, filterXY = TRUE, QCimages = F, plotBMIQ = F)
  # identify MVP and it would be used in DMR calling
  mylimma=champ.MVP(beta.norm = myNorm$beta, pd = data.frame(myLoad$pd), adjPVal = 0.5, adjust.method = "BH",compare.group = c("R","S"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = T)
  # Do DMR calling and save it to result2
  library("doParallel")
  result2=champ.DMR(betaNorm=myNorm$beta,design=myLoad$pd$Sample_Group,maxGap=300,
                    cutoff=0.1,minProbes=3,smooth=TRUE,smoothFunction=loessByCluster,
                    useWeights=FALSE,permutations=NULL,B=10,pickCutoff=FALSE,
                    pickCutoffQ=0.99,nullMethod="bootstrap",verbose=TRUE,cores=3,arraytype="450K",
                    method = "Bumphunter",resultsFile=mylimma,meanLassoRadius=375,minSigProbesLasso=3,minDmrSep=1000,
                    minDmrSize=50,adjPvalProbe=0.5,adjPvalDmr=0.5,pData=myLoad$pd)
  Fileout1<-paste(sampleID,"resultsChamp.tar.gz",sep=".")
  Fileout2<-paste(sampleID,"result.txt",sep=".")
  cmd=paste("tar czvf ",Fileout1," ./resultsChamp",sep="");  
  system(cmd)
#  write.table(result2,file=Fileout2,col.names=NA,row.names=T,sep="\t",quote=F)
  Rlt$dmr=result2
  Rlt$myNorm=myNorm
  Rlt$mylimma=mylimma
  Rlt$myLoad=myLoad
  return(Rlt)
}

# prepare csv samplesheet and each samplesheet will implement the prevous DMR scrpt
drug=read.table("case1.txt",head=T,row.names=1,fill=T,sep="\t")
saminfo=read.table("SampleSheet2016.txt",head=T,fill=T,sep=",",skip=7,as.is=T)
sam2=read.table("Saminfo2.txt",head=T,fill=T,sep="\t")
saminfonew=saminfo[match(sam2[,1],saminfo[,1]),]
saminfonew=data.frame(saminfonew,sam2)

R<-apply(drug,1,function(x) which(x=="R"))
S<-data.frame(apply(drug,1,function(x) which(x=="S")))

for(i in 1:ncol(R)){
  system("rm *csv")
  Res<-match(colnames(drug)[c(R[1,i],R[2,i])],as.character(saminfonew[,ncol(saminfonew)]))
  Sen<-match(colnames(drug)[c(S[1,i],S[2,i])],as.character(saminfonew[,ncol(saminfonew)]))
  saminfonew[Res,3]<-"R"
  saminfonew[Sen,3]<-"S"
  ResSam<-saminfonew[Res,]
  SenSam<-saminfonew[Sen,]
  output=rbind(ResSam,SenSam)
  outputFile<-paste(colnames(R)[i],"samplesheet.csv",sep=".")
  write.table(output,file=outputFile,col.names=T,row.names=F,sep=",",quote=F)
  Rlt<-DMR.CHAMP()
  image=paste(colnames(R)[i],"image.RData",sep=".")
  save.image(file=image)
}
