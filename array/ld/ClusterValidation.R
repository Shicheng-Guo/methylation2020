function (d = NULL, clustering, alt.clustering = NULL, noisecluster = FALSE, 
          silhouette = TRUE, G2 = FALSE, G3 = FALSE, wgap = TRUE, sepindex = TRUE, 
          sepprob = 0.1, sepwithnoise = TRUE, compareonly = FALSE, 
          aggregateonly = FALSE) 
{
  if (!is.null(d)) 
    d <- as.dist(d)
  cn <- max(clustering)
  clusteringf <- as.factor(clustering)
  clusteringl <- levels(clusteringf)
  cnn <- length(clusteringl)
  if (cn != cnn) {
    warning("clustering renumbered because maximum != number of clusters")
    for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
    cn <- cnn
  }
  n <- length(clustering)
  noisen <- 0
  cwn <- cn
  if (noisecluster) {
    noisen <- sum(clustering == cn)
    cwn <- cn - 1
  }
  diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <- within.dist <- between.dist <- numeric(0)
  for (i in 1:cn) cluster.size[i] <- sum(clustering == i)
  pk1 <- cluster.size/n
  pk10 <- pk1[pk1 > 0]
  h1 <- -sum(pk10 * log(pk10))
  corrected.rand <- vi <- NULL
  if (!is.null(alt.clustering)) {
    choose2 <- function(v) {
      out <- numeric(0)
      for (i in 1:length(v)) out[i] <- ifelse(v[i] >= 2, 
                                              choose(v[i], 2), 0)
      out
    }
    cn2 <- max(alt.clustering)
    clusteringf <- as.factor(alt.clustering)
    clusteringl <- levels(clusteringf)
    cnn2 <- length(clusteringl)
    if (cn2 != cnn2) {
      warning("alt.clustering renumbered because maximum != number of clusters")
      for (i in 1:cnn2) alt.clustering[clusteringf == clusteringl[i]] <- i
      cn2 <- cnn2
    }
    nij <- table(clustering, alt.clustering)
    dsum <- sum(choose2(nij))
    cs2 <- numeric(0)
    for (i in 1:cn2) cs2[i] <- sum(alt.clustering == i)
    sum1 <- sum(choose2(cluster.size))
    sum2 <- sum(choose2(cs2))
    pk2 <- cs2/n
    pk12 <- nij/n
    corrected.rand <- (dsum - sum1 * sum2/choose2(n))/((sum1 + 
                                                          sum2)/2 - sum1 * sum2/choose2(n))
    pk20 <- pk2[pk2 > 0]
    h2 <- -sum(pk20 * log(pk20))
    icc <- 0
    for (i in 1:cn) for (j in 1:cn2) if (pk12[i, j] > 0) 
      icc <- icc + pk12[i, j] * log(pk12[i, j]/(pk1[i] * 
                                                  pk2[j]))
    vi <- h1 + h2 - 2 * icc
  }
  if (compareonly) {
    out <- list(corrected.rand = corrected.rand, vi = vi)
  }
  else {
    if (silhouette) 
      require(cluster)
    dmat <- as.matrix(d)
    within.cluster.ss <- 0
    overall.ss <- nonnoise.ss <- sum(d^2)/n
    if (noisecluster) 
      nonnoise.ss <- sum(as.dist(dmat[clustering <= cwn, 
                                      clustering <= cwn])^2)/sum(clustering <= cwn)
    ave.between.matrix <- separation.matrix <- matrix(0, 
                                                      ncol = cn, nrow = cn)
    di <- list()
    for (i in 1:cn) {
      cluster.size[i] <- sum(clustering == i)
      di <- as.dist(dmat[clustering == i, clustering == 
                           i])
      if (i <= cwn) {
        within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
        within.dist <- c(within.dist, di)
      }
      if (length(di) > 0) 
        diameter[i] <- max(di)
      else diameter[i] <- NA
      average.distance[i] <- mean(di)
      median.distance[i] <- median(di)
      bv <- numeric(0)
      for (j in 1:cn) {
        if (j != i) {
          sij <- dmat[clustering == i, clustering == 
                        j]
          bv <- c(bv, sij)
          if (i < j) {
            separation.matrix[i, j] <- separation.matrix[j, 
                                                         i] <- min(sij)
            ave.between.matrix[i, j] <- ave.between.matrix[j, 
                                                           i] <- mean(sij)
            if (i <= cwn & j <= cwn) 
              between.dist <- c(between.dist, sij)
          }
        }
      }
      separation[i] <- min(bv)
      average.toother[i] <- mean(bv)
    }
    average.between <- mean(between.dist)
    average.within <- mean(within.dist)
    nwithin <- length(within.dist)
    nbetween <- length(between.dist)
    between.cluster.ss <- nonnoise.ss - within.cluster.ss
    ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss * 
                                                     (cwn - 1))
    clus.avg.widths <- avg.width <- NULL
    if (silhouette) {
      sii <- silhouette(clustering, dmatrix = dmat)
      sc <- summary(sii)
      clus.avg.widths <- sc$clus.avg.widths
      if (noisecluster) 
        avg.width <- mean(sii[clustering <= cwn, 3])
      else avg.width <- sc$avg.width
    }
    g2 <- g3 <- cn2 <- cwidegap <- widestgap <- sindex <- NULL
    if (G2) {
      splus <- sminus <- 0
      for (i in 1:nwithin) {
        splus <- splus + sum(within.dist[i] < between.dist)
        sminus <- sminus + sum(within.dist[i] > between.dist)
      }
      g2 <- (splus - sminus)/(splus + sminus)
    }
    if (G3) {
      sdist <- sort(c(within.dist, between.dist))
      sr <- nwithin + nbetween
      dmin <- sum(sdist[1:nwithin])
      dmax <- sum(sdist[(sr - nwithin + 1):sr])
      g3 <- (sum(within.dist) - dmin)/(dmax - dmin)
    }
    pearsongamma <- cor(c(within.dist, between.dist), c(rep(0, 
                                                            nwithin), rep(1, nbetween)))
    dunn <- min(separation[1:cwn])/max(diameter[1:cwn], na.rm = TRUE)
    acwn <- ave.between.matrix[1:cwn, 1:cwn]
    dunn2 <- min(acwn[upper.tri(acwn)])/max(average.distance[1:cwn])
    if (wgap) {
      cwidegap <- rep(0, cwn)
      for (i in 1:cwn) if (sum(clustering == i) > 1) 
        cwidegap[i] <- max(hclust(as.dist(dmat[clustering == 
                                                 i, clustering == i]), method = "single")$height)
      widestgap <- max(cwidegap)
    }
    if (sepindex) {
      psep <- rep(NA, n)
      if (sepwithnoise | !noisecluster) {
        for (i in 1:n) psep[i] <- min(dmat[i, clustering != 
                                             clustering[i]])
        minsep <- floor(n * sepprob)
      }
      else {
        dmatnn <- dmat[clustering <= cwn, clustering <= 
                         cwn]
        clusteringnn <- clustering[clustering <= cwn]
        for (i in 1:(n - noisen)) psep[i] <- min(dmatnn[i, 
                                                        clusteringnn != clusteringnn[i]])
        minsep <- floor((n - noisen) * sepprob)
      }
      sindex <- mean(sort(psep)[1:minsep])
    }
    if (!aggregateonly) 
      out <- list(n = n, cluster.number = cn, cluster.size = cluster.size, 
                  min.cluster.size = min(cluster.size[1:cwn]), 
                  noisen = noisen, diameter = diameter, average.distance = average.distance, 
                  median.distance = median.distance, separation = separation, 
                  average.toother = average.toother, separation.matrix = separation.matrix, 
                  ave.between.matrix = ave.between.matrix, average.between = average.between, 
                  average.within = average.within, n.between = nbetween, 
                  n.within = nwithin, max.diameter = max(diameter[1:cwn], 
                                                         na.rm = TRUE), min.separation = sepwithnoise * 
                    min(separation) + (!sepwithnoise) * min(separation[1:cwn]), 
                  within.cluster.ss = within.cluster.ss, clus.avg.silwidths = clus.avg.widths, 
                  avg.silwidth = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma, 
                  dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within/average.between, 
                  ch = ch, cwidegap = cwidegap, widestgap = widestgap, 
                  sindex = sindex, corrected.rand = corrected.rand, 
                  vi = vi)
    else out <- list(n = n, cluster.number = cn, min.cluster.size = min(cluster.size[1:cwn]), 
                     noisen = noisen, average.between = average.between, 
                     average.within = average.within, max.diameter = max(diameter[1:cwn], 
                                                                         na.rm = TRUE), min.separation = sepwithnoise * 
                       min(separation) + (!sepwithnoise) * min(separation[1:cwn]), 
                     ave.within.cluster.ss = within.cluster.ss/(n - noisen), 
                     avg.silwidth = avg.width, g2 = g2, g3 = g3, pearsongamma = pearsongamma, 
                     dunn = dunn, dunn2 = dunn2, entropy = h1, wb.ratio = average.within/average.between, 
                     ch = ch, widestgap = widestgap, sindex = sindex, 
                     corrected.rand = corrected.rand, vi = vi)
  }
  out
}