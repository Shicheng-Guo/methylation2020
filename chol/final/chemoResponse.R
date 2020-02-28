library("TCGAbiolinks")
library("Haplin")
source("GscTools.R")

# 1) receive pid from TCGAbiolinks
pid<-TCGAbiolinks:::getGDCprojects()$project_id
pid<-pid[grep("TCGA",pid)]
# 1) receive pid from github
pid<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/Pid.drugResponse.txt",head=F,sep="\t")
pid<-as.character(pid[,1])

drug2csv<-function(clinical.drug){
  bcr_patient_barcode<-clinical.drug$bcr_patient_barcode
  therapy_types<-clinical.drug$therapy_types
  drug_name<-clinical.drug$drug_name
  measure_of_response<-clinical.drug$measure_of_response
  days_to_drug_therapy_start<-clinical.drug$days_to_drug_therapy_start
  days_to_drug_therapy_end<-clinical.drug$days_to_drug_therapy_end
  therapy_ongoing<-clinical.drug$therapy_ongoing
  new.clinical.drug<-data.frame(bcr_patient_barcode,therapy_types,drug_name,measure_of_response,days_to_drug_therapy_start,days_to_drug_therapy_end,therapy_ongoing)
  return(new.clinical.drug)
}

rlt<-c()
for(i in pid){
  query <- GDCquery(project=i,data.category = "Clinical",file.type = "xml")
  GDCdownload(query)
  clinical.drug <- GDCprepare_clinic(query,"drug")
  drugResponse<-drug2csv(clinical.drug)
  rlt<-rbind(rlt,drugResponse)
}
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/pan")
write.table(rlt,file="Pancancer.drugresponse_20190911.txt",sep="\t",col.names=NA,row.names=T,quote=F)

rlt<-read.table(file="Pancancer.drugresponse_20190911.txt",sep="\t",fill=T,head=T)
rlt<-rlt[,2:ncol(rlt)]
rlt<-subset(rlt,measure_of_response !="" & therapy_types=="Chemotherapy")
rlt<-rlt[match(unique(rlt[,1]),rlt[,1]),]
dim(rlt)
cid<-paste(rlt[,1],"-01",sep="")
methid<-read.table("gdc_manifest.2019-09-11_meth450.txt",head=T)

F1<-id2phen4(methid$filename)[id2phen4(methid$filename) %in% cid]
length(F1)

manifest2barcode("gdc_manifest.2019-09-11_miRNA.txt")  # run in Linux, not works in Windows
system("mv barcode.txt gdc_manifest.2019-09-11_miRNA.barcode.txt")
manifest2barcode("gdc_manifest.2019-09-11_mRNA.txt")  # run in Linux, not works in Windows
system("mv barcode.txt gdc_manifest.2019-09-11_mRNA.barcode.txt")

miRnaid<-read.table("gdc_manifest.2019-09-11_miRNA.barcode.txt",sep="\t",head=T)
F2<-id2phen4(miRnaid$cases.0.samples.0.submitter_id)[id2phen4(miRnaid$cases.0.samples.0.submitter_id) %in% F1]
length(F2)

mRnaid<-read.table("gdc_manifest.2019-09-11_mRNA.barcode.txt",sep="\t",head=T)
F3<-id2phen4(mRnaid$cases.0.samples.0.submitter_id)[id2phen4(mRnaid$cases.0.samples.0.submitter_id) %in% F2]
length(F3)







