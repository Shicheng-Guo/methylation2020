

source("http://www.bioconductor.org/biocLite.R")
biocLite("charm")
biocLite("chipchip")
biocLite("DMRcate")

## script to analyze mm7_ENCODE data
library("charm")
library("DMRcate")
library("oligoClasses")
library("chipchip")
library("oligo")
names=base:::names

####### read in data
## placenta1
dat.pl <- read.xysfiles(c("66138_532.xys","66138_635.xys"))

####### normalization
load("../target.rda")
dat.pl <- quantile.norm(dat.pl,target)

## find blocks
winsize=500 ## 500 probes smoothing
minlen=10000 ## 10kb minimum size 
blk.pl <- blockFinder(dat.pl, span=winsize, lenmin=minlen)

## block list
blklist <- summary(blk.pl)

## pct of regions being blocks
sum(blklist$length)/sum(array.coverage(dat.pl))

## look at blocks
plot(blk.pl)



