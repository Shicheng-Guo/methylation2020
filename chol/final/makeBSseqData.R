makeBSseqData<-function (dat, sampleNames){
    n0 <- length(dat)
    if (missing(sampleNames)) 
        sampleNames <- paste("sample", 1:n0, sep = "")
    alldat <- dat[[1]]
    if (any(alldat[, "N"] < alldat[, "X"], na.rm = TRUE)) 
        stop("Some methylation counts are greater than coverage.\n")
    ix.X <- which(colnames(alldat) == "X")
    ix.N <- which(colnames(alldat) == "N")
    colnames(alldat)[ix.X] <- "X1"
    colnames(alldat)[ix.N] <- "N1"
    if (n0 > 1) {
        for (i in 2:n0) {
            thisdat <- dat[[i]]
            if (any(thisdat[, "N"] < thisdat[, "X"], na.rm = TRUE)) 
                stop("Some methylation counts are greater than coverage.\n")
            ix.X <- which(colnames(thisdat) == "X")
            ix.N <- which(colnames(thisdat) == "N")
            colnames(thisdat)[c(ix.X, ix.N)] <- paste(c("X", 
                "N"), i, sep = "")
            alldat <- merge(alldat, thisdat, all = TRUE)
        }
    }
    ix.X <- grep("X", colnames(alldat))
    ix.N <- grep("N", colnames(alldat))
    alldat[is.na(alldat)] <- 0
    M <- as.matrix(alldat[, ix.X, drop = FALSE])
    Cov <- as.matrix(alldat[, ix.N, drop = FALSE])
    colnames(M) <- colnames(Cov) <- sampleNames
    result <- BSseq(chr = alldat$chr, pos = alldat$pos, M = M, 
        Cov = Cov)
    result
}