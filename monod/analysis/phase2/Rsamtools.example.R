library("Rsamtools")
library("GenomicRanges")
library("GenomicAlignments")

# index CpG location
file<-read.table("/home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt",head=F,sep="\t")
cpg<-list()
for(i in paste("chr",c(1:22,"X","Y","M"),sep="")){
  cpg[[i]]<-file[file[,1] %in% i,2]
  print(i)
}

setwd("/home/shg047/bam")
# quickBamFlagSummary("Indx16_S12.sorted.clipped.bam", main.groups.only=TRUE)
flag1 <- scanBamFlag(isFirstMateRead=NA, isSecondMateRead=NA,isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
param1<-ScanBamParam(flag=flag1, what=c("rname","pos","qwidth","strand","cigar","seq"))
gal1 <- readGAlignments("Indx16_S12.sorted.clipped.bam", use.names=TRUE, param=param1)
data<-mcols(gal1)

data<-mcols(gal1)
strand<-mcols(gal1)$strand
chr<-mcols(gal1)$rname
start<-start(gal1)
end<-end(gal1)
seq<-mcols(gal1)$seq

haplotype<-c()
for(i in 1:length(seq)){
  haplotype.tmp<-haplo(seq[i],chr[i],start[i],end[i],strand[i])
  haplotype<-c(haplotype,haplotype.tmp)
  print(haplotype.tmp)
}
  
haplo<-function(seq,chr,start,end,strand){
  pos<-findcpg(cpg[[chr]],start,end)
  if(length(pos)<2){next}
  offset=pos-start
  if(strand=="+"){
    hap<-paste(unlist(strsplit(as.character(seq),""))[pos],collapse="")
  }else if(strand=="-"){
    hap<-paste(unlist(strsplit(as.character(seq),""))[pos+1],collapse="")
  }
  return(hap)
}

i=5000
myseq=seq[i]
pos=findcpg(cpg[[chr[i]]],start[i],end[i])
hap<-paste(unlist(strsplit(as.character(myseq),""))[pos],collapse="")
hap

x<-c(1:20)
min=10
max=15
findcpg(x,10,15)
findcpg<-function(x,min,max){
i=min
start <-min
end<-max
while(i < max){
  i=i+1
  if(x[i]>=min && x[i-1]<=min){
    start<-i
  }
  if(x[i]<=max){
    end<-i
  }
  if(x[i]>=max && x[i-1]<max){
    i=length(x)
  }
}
return(c(start,end))
}




x<-c(1:100)
findcpg(x,10,20)
is_on_minus <- as.logical(strand(gal1) == "-")
seq[is_on_minus] <- reverseComplement(seq[is_on_minus])
bf <- BamFile("Indx16_S12.sorted.clipped.bam", asMates=F)
paln <- readGAlignmentsList(bf)
j <- junctions(paln, with.revmap=TRUE)
roi <- GRanges("chr1", IRanges(10056, width=500)) 
findOverlaps(paln,roi)
j_overlap <- paln[paln %over% roi]
paln[j_overlap$revmap[[1]]]
? scanBamWhat
names(gal1)
? ScanBamParam
ScanBamParam(which=which, what=what).
? readGAlignmentsList
cigar<-(sort(table(cigar(gal1)), decreasing=TRUE))
cigar[grep("*S*",cigar)]
colSums(cigarOpTable(cigar(gal1)))
table(njunc(gal1))
biocLite("BiSeq")
library(BiSeq)
? BamFile
library("Bismark")
methylation_extractor("Indx16_S12.sorted.clipped.bam")
readBismark()
? readBismark







