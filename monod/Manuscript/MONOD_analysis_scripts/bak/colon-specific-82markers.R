load("./MHL4.RData")
colnames(data)<-gsub(".sorted.clipped.bam.hapInfo.txt.hap","",colnames(data))
colnames(data)<-gsub(".hapInfo.txt.hap","",colnames(data))
saminfo<-read.table("/media/NAS1/shg047/monod/hapinfo/N37Salk.saminfo",sep="\t")
Data=data[,grep("STL|N37|age|ENC|new|centenarian|CTT|HCT|X7.T|X6.T|X6.P|RRBS.6P|X7.P|RRBS.7P|NC.P",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
Data<-rename(Data)
colnames(Data)[grep("age|new|centenarian|WB|middle",colnames(Data))]<-"WBC"
colnames(Data)<-unlist(lapply(colnames(Data),function(x) unlist(strsplit(x,"[.|-]"))[1]))

colon<-bio[which(bio$group=="Colon"),]
input<-Data[match(rownames(colon),rownames(Data)),]
apply(input,2,function(x) mean(x,na.rm=T))
