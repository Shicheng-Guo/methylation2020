

library("DeconRNASeq")
signatures<-matrix(abs(rnorm(121,0.1,0.1)),11,11)
diag(signatures)<-abs(rnorm(11,0.1,0.1))
colnames(signatures)=c("Brain","CCT","Colon","Esophagus","Heart","Intestine","Kidney","Liver","Lung","Stomach","WBC")
VirtualMatrix<-c()
for(F1 in seq(0,0.5,by=0.05)){
  VirtualSample<-F1/2*(signatures[,grep("CCT",colnames(signatures))])+F1/2*(signatures[,grep("Colon",colnames(signatures))])+(1-F1)*(signatures[,grep("WB",colnames(signatures))])
  VirtualMatrix<-cbind(VirtualMatrix,VirtualSample)
}
VirtualMatrix<-data.frame(VirtualMatrix)
colnames(VirtualMatrix)=paste("VM",seq(0,0.5,by=0.01),sep="")

DeconData<-data.frame(VirtualMatrix,signatures)
VirtualMatrix=data.frame(DeconData[,grep("VM",colnames(DeconData))])
Signatures=data.frame(DeconData[,-grep("VM",colnames(DeconData))])
Rlt<-DeconRNASeq(VirtualMatrix,Signatures, checksig=FALSE,known.prop = F, use.scale = TRUE, fig = TRUE)
rlt<-Rlt$out.all
rlt


signatures<-matrix(abs(rnorm(121,0.1,0.1)),11,11)
diag(signatures)<-abs(rnorm(11,0.1,0.1))
colnames(signatures)=c("Brain","CCT","Colon","Esophagus","Heart","Intestine","Kidney","Liver","Lung","Stomach","WBC")
VirtualMatrix<-c()
for(F1 in seq(0,1,by=0.05)){
  VirtualSample<-F1*(signatures[,grep("CCT",colnames(signatures))])+(1-F1)*(signatures[,grep("WB",colnames(signatures))])
  VirtualMatrix<-cbind(VirtualMatrix,VirtualSample)
}
VirtualMatrix<-data.frame(VirtualMatrix)
colnames(VirtualMatrix)=paste("VM",seq(0,1,by=0.05),sep="")

DeconData<-data.frame(VirtualMatrix,signatures)
VirtualMatrix=data.frame(DeconData[,grep("VM",colnames(DeconData))])
Signatures=data.frame(DeconData[,-grep("VM",colnames(DeconData))])
Rlt<-DeconRNASeq(VirtualMatrix,Signatures, checksig=FALSE,known.prop = F, use.scale = TRUE, fig = TRUE)
rlt<-Rlt$out.all
rlt


