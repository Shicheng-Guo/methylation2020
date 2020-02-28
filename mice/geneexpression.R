

for i in `ls *.txt.shore.genes.bed`
do
awk '{print $11,":",$12,"\t",$9,"\t",$14}' OFS="" $i > $i.bin
done 


file=list.files(pattern="*bin")
cor<-c()
for(j in 1:length(file)){
tmp<-read.table(file[j],sep="\t")
cor<-c(cor,as.character(tmp[,1]))  
}

file=list.files(pattern="*bin")
dat<-c()
for(j in 1:length(file)){
  tmp<-read.table(file[j],sep="\t")
  tmp<-tmp[match(names(which(table(cor)==11)),tmp[,1]),]
  dat<-cbind(dat,tmp[,3])  
}

dat<-data.frame(tmp[,2],dat)
colname<-c("gene","HDF" , "hESO7" , "hESO8" , "iPS-S1"  , "iPS-S2"  , "iPS-R1", "iPS-R2", "NT1" ,"NT2", "NT3" ,"NT4")
colnames(dat)<-colname
exp<-read.table("gene_expression.txt",head=T,sep="\t",check.names = F)
exp[2,]
newexp<-exp[,sapply(colnames(dat)[2:12],function(x) match(x,colnames(exp))[1])]

apply(dat[8:26,2:12],1,function(x) cor.test(as.numeric(newexp[2,]),as.numeric(x))$p.value)



