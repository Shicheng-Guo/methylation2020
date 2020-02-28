setwd("C:\\Users\\shg047\\Dropbox\\Project\\methylation\\GastrointestinalCancerProject")
setwd("/home/sguo/Dropbox/Project/methylation/GastrointestinalCancerProject")
files<-list.files()
files
chol<-read.table( "TCGA-CHOL.DMS.txt",head=T)
coad<-read.table("TCGA-COAD.READ.DMS.txt",head=T)
lihc<-read.table("TCGA-LIHC.DMS.txt",head=T)
paad<-read.table("TCGA-PAAD.DMS.txt",head=T)
stad<-read.table("TCGA-STAD.DMS.txt",head=T)

x1<-as.character(chol$V4)
x2<-as.character(coad$V4)
x3<-as.character(lihc$V4)
x4<-as.character(paad$V4)
x5<-as.character(stad$V4)

xx<-c(x1,x2,x3,x4,x5)
yy<-table(xx)[(table(xx))==5]
yy
length(yy)
length(x1)
length(x2)
length(x3)
length(x4)
length(x5)


/media/Home_Raid1/shg047/NAS1/Minghua2017/Gastrointestinal


