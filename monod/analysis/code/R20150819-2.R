v<-c()
for(i in 0:1){
  for(j in 0:1){
    for(k in 0:1){
      for(l in 0:1){
        v<-rbind(v,c(i,j,k,l))
      } 
    } 
  } 
}
v

dim(v)


library("ggplot2")
d <- data.frame(CpG=rep(1:10,4), read=rep(1:4, each=10), value=sample(c(1,16), 40, replace=T))
p<-ggplot() + scale_shape_identity() + geom_point(data=d, aes(x=CpG, y=read, shape=value), size=9)+coord_cartesian(ylim = c(0, 5))
p+geom_line(abline(h=1:4),colour=2)
p
