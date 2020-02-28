BiocManager::install("MEDIPS")
BiocManager::install("gtools")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

BiocManager::install("bigmemory")
BiocManager::install("gplots")

library(MEDIPS)
library(gtools)
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome.Hsapiens.UCSC.hg19)
library(bigmemory)
library(gplots)

#setwd("/gpfs/home/wallen/experiment/experiment/stavros_data")
#setwd("~/Documents/Neuroscience/barnea_lab/rna_seq_experiment/stavros_chip")

setwd("~/Documents/research/methylation/illumina/seqdata")

source("~/r/methyl/LEA/medips.R")
source("~/r/methyl/LEA/util.R")
source("~/r/methyl/LEA/dipdata.R")
source("~/r/seqdata.R")

diffEnrichedROILoRes <- function() {
  pairs <- list(c("omp_medip", "icam_medip"),
                c("ngn_medip", "icam_medip"),
                c("omp_medip", "ngn_medip"),
                c("ngn_hmedip", "icam_hmedip"),
                c("omp_hmedip", "ngn_hmedip"),
                c("omp_hmedip", "icam_hmedip"))
  diffEnrichedFeaturesROILoRes(pairs, "roi/kg_tss.txt", "tss")
  diffEnrichedFeaturesROILoRes(pairs, "roi/kg_genes.txt", "genes")
  diffEnrichedFeaturesROILoRes(pairs, "roi/kg_promoters_10kb.txt", "prom")
  diffEnrichedFeaturesROILoRes(pairs, "roi/kg_exons.txt", "exons")
  diffEnrichedFeaturesROILoRes(pairs, "roi/kg_introns.txt", "introns")
  diffEnrichedFeaturesROILoRes(pairs, "roi/cpgislands.txt", "cpgislands")
}

diffEnrichedFeaturesROILoRes <- function(pairs, roi.path, type) {
  if (!file.exists(roi.path)) {
    quit(paste("roi file at", roi.path, "does not exist", sep=""))
  }
  feat.roi <- loadROI(roi.path)
  colnames(feat.roi) <- c('chr', 'start', 'end')
  path.ex <- paste("_", type, "_diff_lores.txt", sep="")
  foreach (pair = pairs) %dopar% {
    cat(pair[1], pair[2], "\n")
    basename <- paste(pair[1], "_", pair[2], path.ex, sep="") 
    enrich1 <- rawCountsByROILoRes.DipData(load.DipData(pair[1]), feat.roi) 
    enrich2 <- rawCountsByROILoRes.DipData(load.DipData(pair[2]), feat.roi)
    fc <- foldchange(enrich1, enrich2)
    p.val <- diffEnrichmentROI(enrich1, enrich2)
    saveSignificantROI(feat.roi, fc, p.val, paste("diff_enrich_lores/", basename, sep=""), 0.01)
  }  
}

.doROIDiff <- function(pairs, roi, out.rootname="signif") {
  #    foreach (pair = pairs)  %do% {
  for (pair in pairs) {
    fname <- paste(paste("enriched25kb", pair[1], "vs", pair[2], sep="_"), "txt", sep=".")
    cat("loading", pair[1], '\n')
    dd1 <- subsetByROI.DipData(load.DipData(pair[1]), roi)
    cat("loading", pair[2], '\n')
    dd2 <- subsetByROI.DipData(load.DipData(pair[2]), roi)
    cat('computing and saving diff enrich', fname, '\n')
    p.val <- computeAndSaveDiffEnrich.DipData(dd1, dd2, fname)
    cat(sum(p.val < .1), 'locs at less than .1 pval \n')
    if (sum(p.val < .1) > 0) {
      set <- subsetByPos.DipData(dd1, which(p.val < .1))
      writeSubsetROI(mergeROIWindows(set), 1000, paste(paste("roi", out.rootname, sep="/"), pair[1], pair[2], ".txt", sep="_"))
    }
    cat("==============\n")
  }
}

diffEnrichedROI <- function() {
  me.roi <- loadROI("roi/medip_enrich_roi.txt")
  hme.roi <- loadROI("roi/hmedip_enrich_roi.txt")
  
  me.pairs <- list(c("omp_medip_1kb", "icam_medip_1kb"),
                   c("ngn_medip_1kb", "icam_medip_1kb"),
                   c("omp_medip_1kb", "ngn_medip_1kb"))
  
  hme.pairs <- list(c("ngn_hmedip_1kb", "icam_hmedip_1kb"),
                    c("omp_hmedip_1kb", "ngn_hmedip_1kb"),
                    c("omp_hmedip_1kb", "icam_hmedip_1kb"))
  
  .doROIDiff(me.pairs, me.roi)
  .doROIDiff(hme.pairs, hme.roi)
}

diffEnrichedGenesAndPromoters <- function() {
  me.pairs <- list(c("omp_medip_1kb", "icam_medip_1kb"),
                   c("ngn_medip_1kb", "icam_medip_1kb"),
                   c("omp_medip_1kb", "ngn_medip_1kb"))
  hme.pairs <- list(c("ngn_hmedip_1kb", "icam_hmedip_1kb"),
                    c("omp_hmedip_1kb", "ngn_hmedip_1kb"),
                    c("omp_hmedip_1kb", "icam_hmedip_1kb"))
  
  diffEnrichedPromotersROI(hme.pairs, 1000)
  diffEnrichedGenesROI(hme.pairs, 1000)
}

diffEnrichedEverything <- function() {
  me.pairs <- list(c("omp_medip_1kb", "icam_medip_1kb"),
                   c("ngn_medip_1kb", "icam_medip_1kb"),
                   c("omp_medip_1kb", "ngn_medip_1kb"))
  hme.pairs <- list(c("ngn_hmedip_1kb", "icam_hmedip_1kb"),
                    c("omp_hmedip_1kb", "ngn_hmedip_1kb"),
                    c("omp_hmedip_1kb", "icam_hmedip_1kb"))
  #diffEnrichedTSSROI(pairs, 1000)
  diffEnrichedPromotersROI(hme.pairs, 1000, 0.0001)
  diffEnrichedGenesROI(hme.pairs, 1000, 0.0001)
  diffEnrichedExonsROI(hme.pairs, 1000, 0.0001)
  #  diffEnrichedIntronsROI(pairs, 1000)
  #  diffEnrichedCpGIslandsROI(pairs, 1000)
}


diffEnrichedTSSROI <- function(pairs, win.size,alpha=0.1) {
  diffEnrichedFeaturesROI(pairs, win.size, "roi/kg_tss.txt", "tss", alpha=alpha)
}

diffEnrichedIntronsROI <- function(pairs, win.size, alpha=0.1) {
  diffEnrichedFeaturesROI(pairs, win.size, "roi/kg_introns.txt", "introns", alpha=alpha)
}

diffEnrichedExonsROI <- function(pairs, win.size, alpha=0.1) {
  diffEnrichedFeaturesROI(pairs, win.size, "roi/kg_exons.txt", "exons", alpha=alpha)
}

diffEnrichedCpGIslandsROI <- function(pairs, win.size, alpha=0.1) {
  diffEnrichedFeaturesROI(pairs, win.size, "roi/cpgislands.txt", "cpgislands", alpha=alpha)
}

diffEnrichedGenesROI <- function(pairs, win.size, alpha=0.1) {
  diffEnrichedFeaturesROI(pairs, win.size, "roi/kg_genes.txt", "gene", alpha=alpha)
}


diffEnrichedPromotersROI <- function(pairs, win.size, alpha=0.1) {
  diffEnrichedFeaturesROI(pairs, win.size, "roi/kg_promoters_10kb.txt", "prom", alpha=alpha)
}


diffEnrichedFeaturesROI <- function(pairs, win.size, roi.path, type, alpha=0.001) {
  if (!file.exists(roi.path)) {
    quit(paste("roi file at", roi.path, "does not exist", sep=""))
  }
  feat.roi <- loadROI(roi.path)
  path.ex <- paste("_", type, "_diff.txt", sep="")
  diff.out.dir <- paste(type, "_diff/", sep="")
  idx.dd <- load.DipData(pairs[[1]][[1]])
  idx <- computeSubsetIdx.DipData(idx.dd, feat.roi)
  for (pair in pairs) {
    cat(pair[1], pair[2], "\n")
    basename <- paste(pair[1], "_", pair[2], path.ex, sep="")
    fname <- paste(diff.out.dir, basename, sep="")
    dd1 <- subsetByIdx.DipData(load.DipData(pair[1]), idx)
    dd2 <- subsetByIdx.DipData(load.DipData(pair[2]), idx)
    fc <- foldchange(dd1$genome.data[,'raw'], dd2$genome.data[,'raw'])
    p.val <- computeAndSaveDiffEnrich.DipData(dd1, dd2, fname)
    saveSignificantWindowsROI.DipData(dd1, fc, p.val, paste("diff_enrich_roi/", basename, sep=""), window.size=win.size, alpha=alpha)
  }
}

# Writes ROI file of all above threshold enriched windows 
ROIAllEnriched <- function(dsets) {
  all.enriched <- unlist(sapply(dsets, function(i) {
    temp <- load.DipData(i)
    return(whichEnriched.DipData(temp))
  }))
  return(unique(all.enriched))
}

writeEnrichedROI <- function() {
  dsets <- c("omp_medip_25kb", "ngn_medip_25kb", "icam_medip_25kb")
  medip.roi <- ROIAllEnriched(dsets)
  dsets <- c("omp_hmedip_25kb", "ngn_hmedip_25kb", "icam_hmedip_25kb")
  hmedip.roi <- ROIAllEnriched(dsets)
  omp.medip <- load.DipData("omp_medip_25kb")
  writeSubsetROI(omp.medip$genome.data[medip.roi, c('chr', 'pos')], 25000, "roi/medip_enrich_roi.txt")
  writeSubsetROI(omp.medip$genome.data[hmedip.roi, c('chr', 'pos')], 25000, "roi/hmedip_enrich_roi.txt")
}



pairwiseDiffEnrich <- function() {
  #  pairs <- list(c("moe_hmedip_1kb", "moe_ac3_hmedip_1kb"),
  #                c("omp_hmedip_1kb", "icam_hmedip_1kb"),
  #                c("moe_medip_1kb", "moe_ac3_medip_1kb"),
  #                c("omp_medip_1kb", "icam_medip_1kb"))
  pairs <- list(c("omp_hmedip_1kb", "ngn_hmedip_1kb"),
                c("ngn_hmedip_1kb", "icam_hmedip_1kb"),
                c("omp_medip_1kb", "ngn_medip_1kb"),
                c("ngn_medip_1kb", "icam_medip_1kb"))
  
  #  i <- c("moe_medip_1kb", "moe_ac3_medip_1kb")
  foreach(i = pairs) %dopar% {
    fname <- paste(paste(i[1], i[2], sep="_vs_"), "txt", sep=".")
    cat("loading", i[1], '\n')
    dd1 <- load.DipData(i[1])
    cat("loading", i[2], '\n')
    dd2 <- load.DipData(i[2])
    cat('computing and saving diff enrich', fname, '\n')
    computeAndSaveDiffEnrich.DipData(dd1, dd2, fname)
  }                
}

loadMoeMedip <- function() {
  name <- "moe_medip"
  curr.dset.100 <- loadMedips(dsetToPath(name), 100)
  curr.dset.1000 <- loadMedips(dsetToPath(name), 1000)
  curr.dset.10000 <- loadMedips(dsetToPath(name), 10000)
  cat("100bp\n")
  saveMedipsForDipData(curr.dset.100, "moe_medip")
  cat("1kb\n")
  saveMedipsForDipData(curr.dset.1000, "moe_medip_1kb")
  cat("10kb\n")
  saveMedipsForDipData(curr.dset.10000, "moe_medip_10kb")
}

convertAllDipsToDipData25kb <- function() {
  names <- c(MEDIP.DATASETS, HMEDIP.DATASETS)
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 25000)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, paste(name, "25kb", sep="_"))
  }
}

convertAllDipsToDipData100kb <- function() {
  names <- c(MEDIP.DATASETS, HMEDIP.DATASETS)
  foreach(name = names) %dopar% {
    cat(name, '\n')
    curr.dset <- loadMedips(dsetToPath(name), 100000)
    cat("Saving dipdata", name, "\n")
    saveMedipsForDipData(curr.dset, paste(name, "100kb", sep="_"))
  }
}

convertAllDipsToDipData10kb <- function() {
  names <- c("omp_medip", "ngn_medip", "icam_medip", "moe_ac3_medip", "omp_hmedip", "ngn_hmedip", "icam_hmedip", "moe_ac3_hmedip", "moe_hmedip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 10000)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, paste(name, "10kb", sep="_"))
  }
}

convertAllDipsToDipData1kb <- function() {
  names <- c("omp_medip", "ngn_medip", "icam_medip", "moe_ac3_medip", "omp_hmedip", "ngn_hmedip", "icam_hmedip", "moe_ac3_hmedip", "moe_hmedip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 1000)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, paste(name, "1kb", sep="_"))
  }
}

# dont forget to remove random chromosomes before running loadMedips
convertMedipsToDipData <- function() {
  names <- c("omp_medip", "ngn_medip", "icam_medip", "moe_ac3_medip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 100)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, name)
  }
}

convertHMedipsToDipData <- function() {
  names <- c("omp_hmedip", "ngn_hmedip", "icam_hmedip", "moe_ac3_hmedip", "moe_hmedip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading hmedips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 100)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, name)
  }
}

rescaleAndDiff <- function(medips1, medips2, window.size=1000) {
  cat("Rescaling medips1\n")
  medips1.res <- reduceMedipsResolution(medips1, window.size)
  cat("Rescaling medips2\n")
  medips2.res <- reduceMedipsResolution(medips2, window.size)
}

loadAndSaveAllHmeDipData <- function() {
  omp.hmedip <- loadMedips(dsetToPath("omp_hmedip"))
  save(omp.hmedip,file= "omp_hmedip.Rdata")
  ngn.hmedip <- loadMedips(dsetToPath("ngn_hmedip"))
  save(ngn.hmedip, file="ngn.hmedip.Rdata")
  icam.hmedip <- loadMedips(dsetToPath("icam_hmedip"))
  save(icam.hmedip, file="icam_hmedip.Rdata")
  
  moe.ac3.hmedip <- loadMedips(dsetToPath("moe_ac3_hmedip"))  
  save(moe.ac3.hmedip, file="moe_ac3_hmedip.Rdata")
  moe.hmedip <- loadMedips(dsetToPath("moe_hmedip"))
  save(moe.hmedip, file="moe_hmedip.Rdata")
}

selectSignificants <- function(medips.diff, fdr.cutoff=0.01) {
  library(multtest)
  q.values <- mt.rawp2adjp(medips.diff$p.value.wilcox, proc="BH")
  
}
# SCRIPT

files <- c("omp_hmedip.bed", "ngn_hmedip.bed", "icam_hmedip.bed", "moe_hmedip.bed", "moe_ac3_hmedip.bed")
chrs <- c(paste("chr", 1:19, sep=""), "chrX", "chrY")


normByRPM <- function()
  
  autoCorr.dipdata <- function(data_type,corr_type,res,...) {
    require(nlme)
    dipdata.names<-paste(paste(rep(c("omp","ngn","icam","moe","moe_ac3"),each=2),c("medip","hmedip"),sep="_"),res,sep="_")
    setwd("acf")
    quartz("",7,7,file=paste("all",data_type,corr_type,res,"acf.pdf",sep="_"),type="pdf")
    par(mfrow=c(5,2))
    par(mar=c(2,2,2,1)+.1,oma=c(2,2,0,0))
    setwd("..")
    for(dipdata in dipdata.names){
      dip <- load.DipData(dipdata)
      res <- (dip$genome.data[2,"pos"]-dip$genome.data[1,"pos"])/1000
      data <- dip$genome.data[,data_type]
      max <- 20
      title <- paste(unlist(strsplit(dipdata,"_")),collapse=" ")
      data.acf <- acf(data, lag.max = max, type=corr_type, main="",axes=F,xlab="",ylab="",...)
      mtext(title,side=3,line=-1)
      axis(1,at=seq(0,max,10), labels=seq(0,max*res,10*res))
      axis(2)
    }
    mtext("kilobases",side=1,outer=T)
    mtext("ACF",side=2,outer=T)
    
    dev.off()
    
  }

autoCorr.pattern <- function(corr_type,...) {
  require(nlme)
  pattern.names<-paste("mm9_W100_seq_count",c("CG","CA","CC","CT"),sep="")
  setwd("acf")
  quartz("",7,7,file=paste("pattern",corr_type,"acf.pdf",sep="_"),type="pdf")
  par(mfrow=c(2,2))
  par(mar=c(2,2,2,1)+.1,oma=c(2,2,0,0))
  setwd("..")
  for(pattern in pattern.names){
    pattern.data <- LoadPatternCounts(pattern)
    res <- (pattern.data[2,2]-pattern.data[1,2])/1000
    count.data <- pattern.data[,3]
    max <- 60
    title <- strsplit(pattern,"_")[[1]][4]
    acf.data <- acf(count.data, lag.max = max, type=corr_type, main="",axes=F,xlab="",ylab="",...)
    mtext(title,side=3,line=-1)
    axis(1,at=seq(0,max,10), labels=seq(0,max*res,10*res))
    axis(2)
  }
  mtext("kilobases",side=1,outer=T)
  mtext("ACF",side=2,outer=T)
  
  dev.off()
  
}

fft.dipdata <- function(dipdata,data_type,res,...){
  #	dipdata.names<-paste(paste(rep(c("omp","ngn","icam","moe","moe_ac3"),each=2),c("medip","hmedip"),sep="_"),res,sep="_")
  setwd("fft")
  quartz("",7,7,file=paste("all",data_type,res,"acf.pdf",sep="_"),type="pdf")
  #	par(mfrow=c(5,2))
  #	par(mar=c(2,2,2,1)+.1,oma=c(2,2,0,0))
  setwd("..")
  #	for(dipdata in dipdata.names){
  dip <- load.DipData(dipdata)
  res <- (dip$genome.data[2,"pos"]-dip$genome.data[1,"pos"])/1000
  data <- dip$genome.data[,data_type]
  max <- 20
  title <- paste(unlist(strsplit(dipdata,"_")),collapse=" ")
  data.fft <- fft(data)
  magnitude <- Mod(data.fft)
  phase <- Arg(data.fft)
  magnitude_firsthalf <- magnitude[1:(length(magnitude)/2)]
  phase_firsthalf <- phase[1:(length(phase)/2)]
  x.axis <- 1:length(magnitude_firsthalf)/length(magnitude)
  plot(x=x.axis,y=magnitude_firsthalf,type="l")
  mtext(title,side=3,line=-1)
  axis(1,at=seq(0,max,10), labels=seq(0,max*res,10*res))
  axis(2)
  #	}
  mtext("kilobases",side=1,outer=T)
  #	mtext("ACF",side=2,outer=T)
  
  dev.off()
  
}

computeStrandScores <- function(cellType = NULL) {
  #	take plus,minus set for hmc, mc
  #	calculate plus/(plus + minus) score
  #	calculate Pearsons correlation for hmc,mc scores:
  #	not summed though: (Xi - mean(X))(Yi - mean(Y))/(sd(X)*sd(Y))
  calculateScore <- function(plus, minus) {
    return(plus/(plus+minus))
  }	
  calculatePearsons <- function(hmc, mc) {
    return(((hmc-mean(hmc)) * (mc-mean(mc))) / (sd(hmc) * sd(mc)));
  }
  
  mod <- c("hmedip","medip")
  strand <- c("plus","minus")
  
  
  samples <- list()
  for (mod.curr in mod) {
    for (strand.curr in strand) {
      data.path <- paste(cellType,"/",mod.curr,"_",strand.curr,"_1000",sep="")
      list.name <- paste(mod.curr,strand.curr,sep=".")
      samples[[list.name]]<-LoadBM(data.path);
    }
  }
  
  strand.scores <- list()
  for (mod.curr in mod) {
    sample.plus <- samples[[paste(mod.curr,"plus",sep=".")]][,4]
    sample.minus <- samples[[paste(mod.curr,"minus",sep=".")]][,4]
    strand.scores[[mod.curr]] <- calculateScore(sample.plus,sample.minus)
  }
  
  #	Pearsons
  strand.Pearsons <- calculatePearsons(strand.scores[[1]],strand.scores[[2]])
  out.matrix <- as.big.matrix(matrix(cbind(samples[[1]][,1:3],strand.Pearsons)),type="double")
  return(out.matrix)
  
  
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param SeqData 
##' @param annotation 
##' @param seqData.type 
##' @param FUN 
##' @return 
##' @author bcolquitt
generateFeatureProfile <- function(SeqData, annotation, seqData.type = "raw", FUN = mean) {
  dd <- load.DipData(SeqData)
  genome <- dd$genome.data
  chrs <- dd$chr.names
  
  #	anno <- LoadBM(annotation, "anno")
  anno <- read.delim(annotation, header = FALSE, col.names = c("chr","start","end","name","group","strand"), colClasses = c("character","integer","integer","character","integer","character"))
  
  genome.split <- bigsplit(genome, 'chr', splitcol = 2, splitret = "list")
  #        anno.split <- bigsplit(anno, 'chr', splitcol = 2, splitret = "list")
  anno.split <- split(anno, anno$chr)
  features.combined <- foreach(curr.genome = genome.split, curr.anno = anno.split, curr.chr = chrs, .combine = "c") %dopar% {
    features <- rep(NA, times = length(curr.genome))
    match.index <- match(curr.genome, curr.anno$start)         
    #this.anno <- anno[anno[,'chr']==curr.chr,]
    features <- curr.anno$group[match.index]
    return(features)
  }
  
  genome.features <- cbind(genome[,seqData.type], features.combined)
  split.by.features <- bigsplit(genome.features, 2, splitcol = 1)
  
  result <- foreach(group = split.by.features, .combine = "c") %dopar% {
    return(do.call(FUN,list(group)))
  }
  
  return(result)
  
  
  
}


