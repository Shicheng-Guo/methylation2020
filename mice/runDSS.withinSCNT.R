library("DSS")
read_methBED <- function(filename)
{
   meth = read.table(file=filename,stringsAsFactors=F)
   x = meth[,c(1,2,5)]
   x = cbind(x,as.integer(meth[,5]*meth[,4]))
   colnames(x) = c('chr','pos','N','X');
   return(x)
   print(filename)
}
############################################################
#Parameters
fdr_cutoff = 0.05

#Input
setwd("/home/sguo/mice/within")
filenames_1 = c('SCNT_P7C.BED.txt.trim','SCNT_P8B.BED.txt.trim')
samplenames_1 = c('SCNT3','SCNT4')
filenames_2 = c('SCNT_NB3.BED.txt.trim','SCNT_B12.BED.txt.trim')	
samplenames_2 = c('SCNT1','SCNT2')
############################################################
#Read input files
input = list()
filenames=c(filenames_1,filenames_2)
for(file in filenames)
{
  input = append(input,list(read_methBED(file)))
}
save(input,file="input.RData")
#Calling differentially methylated sites
BSobj <- makeBSseqData(input,c("SCNT3","SCNT4","SCNT1", "SCNT2") )
save(BSobj,file="BSobj.RData")
dmlTest <- DMLtest(BSobj, group1=c("SCNT3","SCNT4"), group2=c("SCNT1", "SCNT2"))
save(dmlTest,file="dmlTest.RData")
DMS = callDML(dmlTest)
#Write output to hard disk
fdr_cutoff=0.05
index = which(DMS$fdr <= fdr_cutoff)
write.table(DMS[index,],file="output_DMS_SCNT34-12.tsv",row.names=F,quote=F,sep="\t")

