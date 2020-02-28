gheatmap<-function(datamatrix,rowname=NULL,figuresize=17){
  #usage:
  #row.names= class of the samples such as cancer or normal
  #row: sample
  #column: gene or location
  #install.packages "gplots" if without this package
  
  library("gplots")
    
  color.map <- function(label){ 
    palette(rainbow(length(unique(label))))[array(data.matrix(data.frame(label)))]
  }
  
  if(! is.null(rowname)){
    RowSideColors <- as.character(unlist(lapply(list(rowname), color.map)))
  }else{
    rowname=NULL
  }
  
  data<-as.data.frame(datamatrix)
  #hclustfun = function(x) hclust(x,method = 'centroid')
  pdf("heatmap.pdf",width=figuresize,height=figuresize)  #increase it when fig is big
  heat<-heatmap.2(as.matrix(data),col=topo.colors(75), distfun = function(x) dist(x,method = 'euclidean'),scale="column",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=RowSideColors)  #euclidean
  dev.off()
}

