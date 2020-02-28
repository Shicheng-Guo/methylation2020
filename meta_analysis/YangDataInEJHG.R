

#yang's data: http://www.nature.com/ejhg/journal/v21/n9/full/ejhg2012281a.html

TP2<-c(25,2,48,39,48,86,26,65,64,27,118,48,44,109,131,66,117,46,20)
FP2<-c(4,1,2,2,52,2,1,0,8,2,26,0,1,6,16,1,2,0,26)
FN2<-c(27,10,47,3,63,11,25,13,16,74,0,13,9,61,27,7,62,2,1)
TN2<-c(87,4,36,20,71,50,19,32,18,30,4,11,13,63,14,24,28,29,39)

data2<-data.frame(TP=TP2,FN=FN2,FP=FP2,TN=TN2)
data2
fit.yang<-reitsma(data2)

fit.yang<-reitsma(data2)
fit.guo1<-reitsma(data[input$sampletype=="serum",])
fit.guo2<-reitsma(data[input$sampletype=="tissue",])
summary(fit.guo)
fit.guo<-reitsma(data[input$methods=="MSP",])
summary(fit.guo)
fit.guo<-reitsma(data[input$methods=="qMSP",])
summary(fit.guo)
