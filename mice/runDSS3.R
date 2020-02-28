setwd("/home/sguo/mice/scnt34vips")
library("DSS")
read_methBED <- function(filename)
{
   meth = read.table(file=filename,stringsAsFactors=F,skip=1)
   x = meth[,c(1,2,5)]
   x = cbind(x,as.integer(meth[,5]*meth[,4]))
   colnames(x) = c('chr','pos','N','X');
   return(x)
}
############################################################
#Parameters
fdr_cutoff = 0.05

#Input
filenames_1 = c(
	    'miPS_B3.BED.txt.trim',
	    'miPS_1E12P20.BED.txt.trim','miPS_2A4F1.BED.txt.trim','miPS_2A4F33.BED.txt.trim'
	    )
samplenames_1 = c('miPS2','miPS1E12','miPS2A1','miPS2A33')
filenames_2 = c(
	    'SCNT_P7C.BED.txt.trim',
	    'SCNT_P8B.BED.txt.trim'
	    )		
samplenames_2 = c('SCNT3','SCNT4')
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
BSobj <- makeBSseqData(input,c("miPS2","miPS1E12","miPS2A1","miPS2A33", "SCNT3", "SCNT4") )
save(BSobj,file="BSobj.RData")
dmlTest <- DMLtest(BSobj, group1=c("miPS2", "miPS1E12","miPS2A1","miPS2A33"), group2=c("SCNT3", "SCNT4"))
save(dmlTest,file="dmlTest.RData")
DMS = callDML(dmlTest)
#Write output to hard disk
fdr_cutoff=0.05
index = which(DMS$fdr <= fdr_cutoff)
write.table(DMS[index,],file="output_DMS_4iPS_2SCNT34.tsv",row.names=F,quote=F,sep="\t")

