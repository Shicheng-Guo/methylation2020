grep -P '^cg|Probe\sid\s\(' DNA__Illumina_450K_methylation_Beta_values.txt.txt > DNA__Illumina_450K_methylation_Beta_values.txt
data<-read.table("DNA__Illumina_450K_methylation_Beta_values.txt",sep="\t",head=T,row.names=1)
newdata<-data[,13:(ncol(data)-1)]
write.table(colnames(newdata),file="NCI60-meth450k.cells.txt",sep="\t",quote=F)
output=unique(newdata)
RowName<-unlist(lapply(rownames(output),function(x) unlist(strsplit(x,"~"))[1]))
rownames(output)=RowName
nci60methdata=output
save(nci60methdata,file="NCI60-meth450k.beta.RData")
write.table(output,file="NCI60-meth450k.beta.txt",sep="\t",quote=F,col.names = NA,row.names = T)

cellLines<-names(table(unlist(lapply(colnames(nci60methdata),function(x) unlist(strsplit(x,"[.]"))[1]))))
for(i in cellLines){
  beta<-nci60methdata[,grep(i,colnames(nci60methdata))]
  save(beta,file=paste(i,"beta.RData",sep="."))
  print(i)
}




grep -P '^cg|Probe\sid\s\(' DNA__Illumina_450K_methylation_Beta_values.txt.txt > DNA__Illumina_450K_methylation_Beta_values.txt
data<-read.table("DNA__Illumina_450K_methylation_Beta_values.txt",sep="\t",head=T,row.names=1)

CHR<-paste("chr",data[,3],sep="")
START<-data[,4]-1
END<-data[,4]
for(i in 13:(ncol(data)-1)){
MValue<-round(as.numeric(as.character(data[,i]))*100)
UValue<-round((1-as.numeric(as.character(data[,i])))*100)
rlt<-data.frame(CHR,START,END,MValue,UValue)
filename<-gsub("[.]","-",colnames(data)[i])
write.table(rlt,file=paste(filename,".hg19.bed",sep=""),col.names = F,row.names = F,quote = F,sep="\t")
print(i)
}


  
newdata<-data[,13:(ncol(data)-1)]
write.table(colnames(newdata),file="NCI60-meth450k.cells.txt",sep="\t",quote=F)
output=unique(newdata)
RowName<-unlist(lapply(rownames(output),function(x) unlist(strsplit(x,"~"))[1]))
rownames(output)=RowName
nci60methdata=output
save(nci60methdata,file="NCI60-meth450k.beta.RData")
write.table(output,file="NCI60-meth450k.beta.txt",sep="\t",quote=F,col.names = NA,row.names = T)

cellLines<-names(table(unlist(lapply(colnames(nci60methdata),function(x) unlist(strsplit(x,"[.]"))[1]))))
for(i in cellLines){
  beta<-nci60methdata[,grep(i,colnames(nci60methdata))]
  save(beta,file=paste(i,"beta.RData",sep="."))
  print(i)
}


