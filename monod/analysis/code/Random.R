

x<-c(3,4,4,4,5,5,3,6,5,4)
y<-matrix(0,10,10)
y

for(j in 1:500){
x1<-matrix((rnorm(100,0.1,0.3)),10,10)
x1[x1<0.6]<-0
x1[x1>0.6]<-1
for(i in 1:10){
  x1[i,i]<-x[i]-sum(x1[-i,i])
}
y=y+x1
}
y1<-y/500
y1

z<-c()
for(i in 1:10){
  z<-c(z,y1[i,i]/x[i])
}
z

rlt<-rbind(y1,z)
write.table(rlt,file="rlt.txt",sep="\t",col.names=NA,row.names=T,quote=F)

getwd()
x<-c(30,29,10,35,40)
y<-matrix(0,10,5)
y

for(j in 1:500){
  x1<-matrix(abs(rnorm(50,rnorm(1,0.1,0.5),abs(rnorm(1,0.1,0.5)))),10,5)
  x1[x1<0.6]<-0
  x1[x1>0.6]<-1
  for(i in 1:5){
    x1[i,i]<-x[i]-sum(x1[-i,i])
  }
  y=y+x1
}
y1<-y/500
y1

z<-c()
for(i in 1:5){
  z<-c(z,y1[i,i]/x[i])
}
z


rlt<-rbind(y1,z)
write.table(rlt,file="rlt-may242016.txt",sep="\t",col.names=NA,row.names=T,quote=F)
getwd()

