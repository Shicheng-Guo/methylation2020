


./liftOver H1.LMR.hg18.bed hg18ToHg19.over.chain H1.LMR.hg19.bed tmp
./liftOver H1.UMR.hg18.bed hg18ToHg19.over.chain H1.UMR.hg19.bed tmp
./liftOver ME.LMR.hg18.bed hg18ToHg19.over.chain ME.LMR.hg19.bed tmp
./liftOver ME.UMR.hg18.bed hg18ToHg19.over.chain ME.UMR.hg19.bed tmp
./liftOver NPC.LMR.hg18.bed hg18ToHg19.over.chain NPC.LMR.hg19.bed tmp
./liftOver NPC.UMR.hg18.bed hg18ToHg19.over.chain NPC.UMR.hg19.bed tmp
./liftOver MSC.LMR.hg18.bed hg18ToHg19.over.chain MSC.LMR.hg19.bed tmp
./liftOver MSC.UMR.hg18.bed hg18ToHg19.over.chain MSC.UMR.hg19.bed tmp
./liftOver TBL.LMR.hg18.bed hg18ToHg19.over.chain TBL.LMR.hg19.bed tmp
./liftOver TBL.UMR.hg18.bed hg18ToHg19.over.chain TBL.UMR.hg19.bed tmp


./liftOver H1.DMV.hg18.bed hg18ToHg19.over.chain H1.DMV.hg19.bed tmp
./liftOver ME.DMV.hg18.bed hg18ToHg19.over.chain ME.DMV.hg19.bed tmp
./liftOver NPC.DMV.hg18.bed hg18ToHg19.over.chain NPC.DMV.hg19.bed tmp
./liftOver MSC.DMV.hg18.bed hg18ToHg19.over.chain MSC.DMV.hg19.bed tmp
./liftOver TBL.DMV.hg18.bed hg18ToHg19.over.chain TBL.DMV.hg19.bed tmp


Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

bed<-read.table("WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed",sep="\t",as.is=T)
  
x1<-read.table("/home/sguo/annotation/H1.UMR.hg19.bed",sep="\t",as.is=T)
x2<-read.table("/home/sguo/annotation/H1.LMR.hg19.bed",sep="\t",as.is=T)
x3<-read.table("/home/sguo/annotation/ME.UMR.hg19.bed",sep="\t",as.is=T)
x4<-read.table("/home/sguo/annotation/ME.LMR.hg19.bed",sep="\t",as.is=T)
x5<-read.table("/home/sguo/annotation/NPC.UMR.hg19.bed",sep="\t",as.is=T)
x6<-read.table("/home/sguo/annotation/NPC.LMR.hg19.bed",sep="\t",as.is=T)
x7<-read.table("/home/sguo/annotation/MSC.UMR.hg19.bed",sep="\t",as.is=T)
x8<-read.table("/home/sguo/annotation/MSC.LMR.hg19.bed",sep="\t",as.is=T)
x9<-read.table("/home/sguo/annotation/TBL.UMR.hg19.bed",sep="\t",as.is=T)
x10<-read.table("/home/sguo/annotation/TBL.LMR.hg19.bed",sep="\t",as.is=T)

x11<-read.table("/home/sguo/annotation/H1.DMV.hg19.bed",sep="\t",as.is=T)
x12<-read.table("/home/sguo/annotation/ME.DMV.hg19.bed",sep="\t",as.is=T)
x13<-read.table("/home/sguo/annotation/NPC.DMV.hg19.bed",sep="\t",as.is=T)
x14<-read.table("/home/sguo/annotation/MSC.DMV.hg19.bed",sep="\t",as.is=T)
x15<-read.table("/home/sguo/annotation/TBL.DMV.hg19.bed",sep="\t",as.is=T)



y1<-Rbedtools(functionstring="intersectBed",bed1=x1,bed2=bed,opt.string="-wa -u")
y2<-Rbedtools(functionstring="intersectBed",bed1=x2,bed2=bed,opt.string="-wa -u")
y3<-Rbedtools(functionstring="intersectBed",bed1=x3,bed2=bed,opt.string="-wa -u")
y4<-Rbedtools(functionstring="intersectBed",bed1=x4,bed2=bed,opt.string="-wa -u")
y5<-Rbedtools(functionstring="intersectBed",bed1=x5,bed2=bed,opt.string="-wa -u")
y6<-Rbedtools(functionstring="intersectBed",bed1=x6,bed2=bed,opt.string="-wa -u")
y7<-Rbedtools(functionstring="intersectBed",bed1=x7,bed2=bed,opt.string="-wa -u")
y8<-Rbedtools(functionstring="intersectBed",bed1=x8,bed2=bed,opt.string="-wa -u")
y9<-Rbedtools(functionstring="intersectBed",bed1=x9,bed2=bed,opt.string="-wa -u")
y10<-Rbedtools(functionstring="intersectBed",bed1=x10,bed2=bed,opt.string="-wa -u")

y11<-Rbedtools(functionstring="intersectBed",bed1=x11,bed2=bed,opt.string="-wa -u")
y12<-Rbedtools(functionstring="intersectBed",bed1=x12,bed2=bed,opt.string="-wa -u")
y13<-Rbedtools(functionstring="intersectBed",bed1=x13,bed2=bed,opt.string="-wa -u")
y14<-Rbedtools(functionstring="intersectBed",bed1=x14,bed2=bed,opt.string="-wa -u")
y15<-Rbedtools(functionstring="intersectBed",bed1=x15,bed2=bed,opt.string="-wa -u")


r1<-nrow(y1)/nrow(x1)
r2<-nrow(y2)/nrow(x2)
r3<-nrow(y3)/nrow(x3)
r4<-nrow(y4)/nrow(x4)
r5<-nrow(y5)/nrow(x5)
r6<-nrow(y6)/nrow(x6)
r7<-nrow(y7)/nrow(x7)
r8<-nrow(y8)/nrow(x8)
r9<-nrow(y9)/nrow(x9)
r10<-nrow(y10)/nrow(x10)

r11<-nrow(y11)/nrow(x11)
r12<-nrow(y12)/nrow(x12)
r13<-nrow(y13)/nrow(x13)
r14<-nrow(y14)/nrow(x14)
r15<-nrow(y15)/nrow(x15)

paste(r11,r12,r13,r14,r15,sep=",")

lda<-read.table("/home/sguo/annotation/LAD.hg19.bed",sep="\t",as.is=T)
locks<-read.table("/home/sguo/annotation/LOCK.hg19.bed",sep="\t",as.is=T)
cpgi<-read.table("/home/sguo/annotation/CpGI.hg19.bed",sep="\t",as.is=T)

lda2<-Rbedtools(functionstring="intersectBed",bed1=lda,bed2=bed,opt.string="-wa -u")
locks2<-Rbedtools(functionstring="intersectBed",bed1=locks,bed2=bed,opt.string="-wa -u")
cpgi2<-Rbedtools(functionstring="intersectBed",bed1=cpgi,bed2=bed,opt.string="-wa -u")

rlda<-nrow(lda2)/nrow(lda)
rlocks<-nrow(locks2)/nrow(locks)
rcpgi<-nrow(cpgi2)/nrow(cpgi)


x11<-read.table("/home/sguo/annotation/H1.DMV.hg19.bed",sep="\t",as.is=T)
x12<-read.table("/home/sguo/annotation/ME.DMV.hg19.bed",sep="\t",as.is=T)
x13<-read.table("/home/sguo/annotation/NPC.DMV.hg19.bed",sep="\t",as.is=T)
x14<-read.table("/home/sguo/annotation/MSC.DMV.hg19.bed",sep="\t",as.is=T)
x15<-read.table("/home/sguo/annotation/TBL.DMV.hg19.bed",sep="\t",as.is=T)


paste(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,sep=",")
x<-c(0.888,0.48,0.27,0.28,0.26,0.26,0.26,0.076,0.046,0.0418,0.0474,0.036,0.054,0.031,0.071,0.04,0.052,0.036)
names(x)<-c("LDA","LOCKS","H1 DMR","ME DMR","NPC DMR","MSC DMR","TBL DMR","CpGI","H1 UMR","H1 LMR","ME UMR","ME LMR","NPC UMR","NPC LMR","MSC UMR","MSC LMR","TBL UMR","TBL LMR")
barplot(x, col =rainbow(12),horiz=F,cex.names = 0.7,ylim=c(0,1))

