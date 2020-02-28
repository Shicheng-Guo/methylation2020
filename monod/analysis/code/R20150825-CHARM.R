source("http://bioconductor.org/biocLite.R")
biocLite("charm")
library("charm")

d1<-read.table("wgEncodeDukeMapabilityRegionsExcludable.bed",sep="\t")
d2<-read.table("wgEncodeDacMapabilityConsensusExcludable.bed",sep="\t")