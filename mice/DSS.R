


data<-read.table("/home/zhl002/ZL_LTS33/mouse_WGBS/bedfiles/Mice.iPS.ES.shore.bed",head=F,sep="\t",as.is=T)
newfdr=p.adjust(data[,4],method="fdr")
newdata<-cbind(data,newfdr)

sig<-subset(newdata,newfdr<0.05)
nonsig<-subset(newdata,newfdr>=0.05)

nrow(sig)
nrow(nonsig)

write.table(sig,file="mice.shore.newfdr.sig.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(nonsig,file="mice.shore.newfdr.nonsig.txt",sep="\t",quote=F,row.names=F,col.names=F)


cp /home/shg047/reprogramming/mice.shore.newfdr.sig.txt ./
cp /home/shg047/reprogramming/mice.shore.newfdr.nonsig.txt ./
  
  
  
cp ../miPS_B3.BED.txt.trim ./
cp ../miPS_1E12P20.BED.txt.trim ./
cp ../miPS_2A4F1.BED.txt.trim ./
cp ../miPS_2A4F33.BED.txt.trim ./
cp ../SCNT_P7C.BED.txt.trim ./
cp ../SCNT_P8B.BED.txt.trim ./


  