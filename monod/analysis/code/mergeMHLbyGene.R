

d<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt.sort.bed.gene.trim.bed",as.is=T)
unigene<-sort(unique(d[,5]))[2:length(unique(d[,5]))]

data1<-c()
for(i in 1:length(unigene)){
nrow<-match(unigene[i],d[,5])  
data1<-rbind(data1,c(d[nrow[1],5],colSums(d[nrow,6:ncol(d)])))
print(i)
}
data1<-data.frame(data1)

nrow<- which(d[,5]==".")
data2<-data.frame(d[nrow,4],d[nrow,6:ncol(d)])

A1<-(matrix(unlist(data1[,2:ncol(data1)]), ncol = ncol(data1)-1, byrow = F))
A2<-(matrix(unlist(data2[,2:ncol(data2)]), ncol = ncol(data2)-1, byrow = F))
rownames(A1)<-data1[,1]
rownames(A2)<-data2[,1]
data<-data.frame(rbind(A1,A2))

# find the colum information
setwd("/home/shg047/monod/mhl")
sam<-read.table("../sampleinfo.sort.txt",sep="\t",as.is=T)
new<-cbind(sapply(sam[,1],function(x) unlist(strsplit(x,":"))[2]),sam[,2])
new<-data.frame(new)
file<-list.files(pattern="*stringent_mhl_matrix.txt")[1:6]
file1<-read.table(file[1],head=T,sep="\t",row.names=1,as.is=T,check.names=F)
cor1<-match(colnames(file1),new[,1])
lab1<-new[cor1,2]
u<-which(lab1=="Cancer")
v<-which(lab1=="Normal")


data1=matrix(as.numeric(as.matrix(data)),ncol=ncol(data),byrow=F)
keep2<-c()
for(i in 1:nrow(data1)){
  print(i)
  p1<-sum(data1[i,u]>0)
  p2<-sum(data1[i,v]==0)
  if(p2==8 & p1>0){
    tmp<-c(i,p1)
    keep2<-rbind(keep2,tmp)
  }  
}
max(keep2[,2])

v1<-apply(file1,1,function(x) sum(x[v]==0))
length(which(v1==8))
u1<-apply(file1[v1,],1,function(x) sum(x[v]==0))

wilcox.test(x[u],x[v])$p.value
p.adjust<-p.adjust(p, "BH") 

pdf("var.pdf")
hist(var,breaks=200)
dev.off()
d<-read.table("1407-combined_RRBS_mld_blocks_stringent_mhl_matrix.txt.bed",sep="\t",as.is=T)









