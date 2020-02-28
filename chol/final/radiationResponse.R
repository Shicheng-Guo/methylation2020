radiation2csv<-function(clinical.radiation){
  bcr_patient_barcode<-clinical.radiation$bcr_patient_barcode
  radiation_type<-clinical.radiation$radiation_type
  radiation_dosage<-clinical.radiation$radiation_dosage
  radiation_treatment_ongoing<-clinical.radiation$radiation_treatment_ongoing
  measure_of_response<-clinical.radiation$measure_of_response
  days_to_radiation_therapy_start<-clinical.radiation$days_to_radiation_therapy_start
  days_to_radiation_therapy_end<-clinical.radiation$days_to_radiation_therapy_end
  measure_of_response<-clinical.radiation$measure_of_response
  new.clinical.radiation<-data.frame(bcr_patient_barcode,radiation_type,radiation_dosage,radiation_treatment_ongoing,days_to_radiation_therapy_start,days_to_radiation_therapy_end,measure_of_response)
  return(new.clinical.radiation)
}

setwd("~/hpc/project/TCGA")
library("TCGAbiolinks")
pid<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/Pid.drugResponse.txt",head=F,sep="\t")
Radiation<-c()
for(i in 1:nrow(pid)){
  query <- GDCquery(project =pid[i,1],data.category = "Clinical", file.type = "xml")
  GDCdownload(query)
  clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
  radiation<-unique(radiation2csv(clinical.radiation))
  Radiation<-rbind(Radiation,radiation)
}
write.table(Radiation,file="Pancancer.radiationresponse.txt",sep="\t",col.names = NA,row.names = T,quote=F)


x<-c(0,0,0,1,1,1)
y<-c(0,0,0,0,1,1)
cor(x,y)

x<-c(1,0,0,1,1,0)
y<-c(1,0,0,0,1,0)
cor(x,y)

gdc_manifest.2019-09-12_image.txt
