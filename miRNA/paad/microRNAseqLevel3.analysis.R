# this is special for paad analysis
setwd("/home/sguo/paad/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3")
file=list.files(pattern="*mirna.quantification.txt")
sampleinfo=read.table("/home/sguo/paad/file_manifest.txt",head=T,sep="\t")
sampleinfo2<-read.table("/home/sguo/paad/Clinical/Biotab/nationwidechildrens.org_clinical_patient_paad.txt",head=T,sep="\t")

sta1<-which(sampleinfo2$clinical_stage=="Stage I")
sta2<-which(sampleinfo2$clinical_stage=="Stage II")
ssa<-paste(sampleinfo2[c(sta1,sta2),1],"-01",sep="")
xx<-as.numeric(sapply(ssa,function(x) grep(x,file)))
xx<-xx[-which(is.na(xx))]
file<-file[xx]
sampleid=sampleinfo[match(file,sampleinfo[,ncol(sampleinfo)]),5]
data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data<-cbind(data,tmp[,2])
  print(i)
}
rownames(data)=tmp[,1]
colnames(data)=sampleid
data<-data.matrix(data)
data2<-data[-which(apply(data,1,function(x){sum(is.na(x))})>0),]
dat<-t(data)
result<-list()
j=0
rat=5000000
for (i in 1:50){
  for(feature in seq(10,300,by=20)){
    print(j)
    j=j+1
    result[[j]]<-Randomized_mcluster_rep(dat,3,feature)  
    real<-sampleinfo2[match(substr(rownames(dat),1,12),sampleinfo2[,1]),]$clinical_stage
    clus<-result[[j]]$cluster
    rlt<-table(real,clus)
    print(rlt)
    ratmp<-(rlt[2,1]*rlt[2,2]*rlt[2,2])+(rlt[3,1]+rlt[3,2]+rlt[3,3])
    if(ratmp<rat){
      rat<-ratmp
      print(i)
      print(rlt)
    }
  }
}






## function from Nan Lin, you need source these function firstly. 
##################################load the require functions#################
##########################randomized feature selection with number of feature input#######
library(irlba)
library(MASS)
Randomized_fsk <- function(A,k=num_cluster,r=num_feature){
  m <- dim(A)[1]
  n <- dim(A)[2]
  R <- matrix(rnorm(n*r,0,1),nrow=n)
  Y <- as.matrix(A) %*% R
  Q <- qr.Q(qr(A))
  
  
  tQA_svd <- irlba(t(Q) %*% as.matrix(A),nu=k,nv=k,adjust=3)
  Z <- tQA_svd$v[,1:k]
  
  ####calcualte the P vector######
  P <- as.vector(rowSums(Z^2)/sum(diag(Z %*% t(Z))))
  
  omega <- matrix(rep(0,n*r),nrow=n)
  S <- matrix(rep(0,r*r),nrow=r)
  
  id <- sample(c(1:n),r, replace=TRUE, prob=P)
  for (i in 1:r){
    omega[id[i],i] <- 1 
    S[i,i] <- 1/sqrt(i*P[id[i]])
  }
  C <- as.matrix(A) %*% omega %*% S
  ###define the objective function to qualify the goodness of fit
  A_svd <- irlba(A,nu=k,nv=k,adjust=3)
  diff_ZV <- sum(diag((A_svd$v - Z) %*% t(A_svd$v - Z)))
  
  return(list("C"= C,"feature"=id,"objective"=diff_ZV,"weight" = diag(S)))
}


Randomized_cluster <- function(A,num_cluster,num_feature){
  tem <- Randomized_fsk(A,num_cluster,num_feature)
  #cluster <- clusterout_01(kmeans(tem[[1]],2)$cluster,status)
  out <- list("data"=tem[[1]],"feature"=tem[[2]],"objective"=tem[[3]])
  return(out)
}

######################multiple cluster by randomization############
Randomized_mcluster_rep <- function(A,num_cluster,num_feature){
  n <- dim(A)[2]
  data <- list()
  feature <- matrix(rep(0,n*num_feature),nrow=n)
  objective <- matrix(rep(0,n),nrow=n)
  acu <- matrix(rep(0,n),nrow=n)
  table <- list()
  
  for (i in 1:n){
    transition <- Randomized_cluster(A,num_cluster,num_feature)    
    data[[i]] <- transition[[1]]
    feature[i,] <- transition[[2]]
    objective[i,] <- transition[[3]]
    #acu[i] <- transition[[4]]
    #obj_acu_fea <- cbind(objective,acu,feature)
    #table[[i]] <- transition[[5]]  
  }
  selection <- matrix(rep(0,n*2),nrow=n)
  for (i in 1:n){
    selection[i,1] <- i
    selection[i,2] <- sum(feature==i)
  }
  sidm <- selection[order(selection[,2],decreasing=TRUE),]
  sid <- sidm[1:num_feature,1]
  cluster <- kmeans(A[,sid],num_cluster)$cluster
  #out <- clusterout_01(kmeans(A[,sid],num_cluster)$cluster,status)
  #out <- list("parameter"=obj_acu_fea,"data"=data,"table"=table)
  return(list("cluster"=cluster,"id"=sid))
}

################################################################
################################################################
###########functions to pick up the best accuracy###############
################################################################
################################################################
permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

opt_class <- function(X,n){
  Permu <- permutations(n)
  acu <- rep(0,factorial(n)) 
  for (i in 1:factorial(n)){
    acu[i] <- sum(diag(X[Permu[i,],]))
  }
  return(max(acu))
}






