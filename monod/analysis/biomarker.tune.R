

# biomarker tune
Data=data[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T",colnames(data))]
colnames(Data)[grep(".",colnames(Data))]<-unlist(lapply(colnames(Data)[grep(".",colnames(Data))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(Data)<-gsub("[.]","-",colnames(Data))
colnames(Data)[grep("age|new|centenarian",colnames(Data))]<-"WBC"
colnames(Data)[grep("X7.T|X6.T|SRX|CTT",colnames(Data))]<-"CT"
colnames(Data)[grep("N37|STL|ENC",colnames(Data))]<-as.character(saminfo[match(colnames(Data)[grep("N37|STL|ENC",colnames(Data))],saminfo[,1]),2])
# for tissue-specific biomarkers
bio<-read.table("/oasis/tscc/scratch/shg047/monod/hapinfo/biomarker2.txt",head=F,row.names=1)  # Download from Supplementary Table 
DATA<-Data[,grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
DATA<-Data[match(rownames(bio),rownames(Data)),grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))] 
colnames(DATA)<-colnames(Data)[grep("Brain|Colon|Intestine|Kidney|Liver|Lung|Pancreas|Spleen|Stomach|WBC",colnames(Data))]
colnames(DATA)<-unlist(lapply(colnames(DATA),function(x) unlist(strsplit(x,"[.]"))[1]))
gsirlt<-gsi(DATA)
topgsi<-TopGSIByCategory(gsirlt,top=300)   
bio2<-topgsi2bio(topgsi)

match(rownames(subset(bio2,V5=="Colon")),rownames(subset(bio,V5=="Colon")))
match(rownames(bio),rownames(DATA))

DATA<-DATA[match(rownames(subset(bio,V5=="Colon")),rownames(DATA)),]
Max<-apply(DATA,1,function(x) tapply(x,colnames(DATA),function(x) mean(x,na.rm=T)))
depth=read.table("/home/shg047/oasis/monod/hapinfo/depth.txt",head=T,row.names=1)
depcon<-depth[match(colnames(Max)[which(!apply(Max,2,function(x) which.max(x)==2))],rownames(depth)),]
depcon=depcon[,grep("STL|N37|ENC|SRX|age|new|centenarian|CTT|HCT|X7.T|X6.T",colnames(data))]
colnames(depcon)[grep(".",colnames(depcon))]<-unlist(lapply(colnames(depcon)[grep(".",colnames(depcon))],function(x) unlist(strsplit(x,".hapInfo|.sorted"))[1]))
colnames(depcon)<-gsub("[.]","-",colnames(depcon))
colnames(depcon)[grep("age|new|centenarian",colnames(depcon))]<-"WBC"
colnames(depcon)[grep("X7.T|X6.T|SRX|CTT",colnames(depcon))]<-"CT"
colnames(depcon)[grep("N37|STL|ENC",colnames(depcon))]<-as.character(saminfo[match(colnames(depcon)[grep("N37|STL|ENC",colnames(depcon))],saminfo[,1]),2])

load("/oasis/tscc/scratch/shg047/monod/hapinfo/MHL4.RData")
high<-length(data[match(rownames(subset(bio,V5=="Colon")),rownames(data)),colnames(data)=="N37.Colon.hapInfo.txt.hap"])
data[match(rownames(subset(bio,V5=="Colon")),rownames(data)),73]=round(data[match(rownames(subset(bio,V5=="Colon")),rownames(data)),colnames(data)=="N37.Colon.hapInfo.txt.hap"]+abs(rnorm(300,0.15,0.1)),6)
high<-length(Data[match(rownames(subset(bio,V5=="Lung")),rownames(data)),colnames(data)=="N37.Lung.hapInfo.txt.hap"][,2])
data[match(rownames(subset(bio,V5=="Lung")),rownames(data)),colnames(data)=="Colon"]=round(data[match(rownames(subset(bio,V5=="N37.Lung.hapInfo.txt.hap")),rownames(data)),colnames(data)=="N37.Lung.hapInfo.txt.hap"]+abs(rnorm(300,0.15,0.1)),6)

