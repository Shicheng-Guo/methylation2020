library("grDevices")
plot(1:20,pch=20,col=col)
col=colorRampPalette(c("white", "red"))(20) 
M <- matrix(runif(100),10,10)
M[lower.tri(M)] <- NA
image(M,col = col,frame=F,xaxt="n",yaxt="n")

library("ggplot2")
d <- data.frame(CpG=rep(1:10,4), read=rep(1:4, each=10), value=sample(c(1,16), 40, replace=T))
ggplot() + scale_shape_identity() + geom_point(data=d, aes(x=CpG, y=read, shape=value), size=9)+coord_cartesian(ylim = c(0, 5))





