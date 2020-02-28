

###########################################################################
# This version R program implementation by Futao Zhang 07/29/2013
###########################################################################
library(MASS)
library(Matrix)
library(fda)
library(mvtnorm)

mvFRG_interaction=function(phenoData,x_A,x_B)
{	
	phenoData=as.matrix(phenoData)
	sample_num = dim(x_A)[1]
	basis_num_A = dim(x_A)[2]
	basis_num_B = dim(x_B)[2]
			
     gamma<- matrix(0,sample_num,basis_num_A*basis_num_B)

     for(i_ in 1:sample_num)
		for(j_ in 1:basis_num_A)
			for(k_ in 1:basis_num_B)
                 	gamma[i_,((j_-1)*basis_num_B)+k_]<-x_A[i_,j_]*x_B[i_,k_]

     W<-cbind(x_A,x_B,gamma)
     WTW=t(W)%*%W
     iWTW=ginv(WTW)
	iWTW_WT=iWTW%*%t(W)
     b_hat=iWTW_WT%*%phenoData
     Y_WB=phenoData-W%*%b_hat
      

 	big_sigma_hat=t(Y_WB)%*%Y_WB/dim(phenoData)[1]
	   
     capital_A= iWTW_WT[-1:-(basis_num_A +basis_num_B),]
	if(is.null(dim(capital_A))) capital_A=t(as.matrix(capital_A))
	
     little_gamma_hat=capital_A%*%phenoData
	A_AT=capital_A%*%t(capital_A)

	big_lamda=kronecker(big_sigma_hat,A_AT)
      
    T=as.numeric(t(as.matrix(as.vector(little_gamma_hat)))%*%ginv(big_lamda)%*%as.matrix(as.vector(little_gamma_hat)))
	
	#rk=qr(big_lamda)$rank
	eigenval<-svd(var(big_lamda))$d
    rk=length(which(eigenval>=1e-8))
  if(rk==0) rk=1
    rlt=pchisq(T,rk,lower.tail=F)
	rlt=c(rlt,rk)
    rlt
}



load("expansion.RData")
load("pair_list.RData")

#names(geno_info)
spl_num=dim(pheno_info)[1]
rng=0

pheno <-  pheno_info[,-1:-3]	
pheno <-  log(pheno)

#c=0.5
#pheno[,1]=qnorm((rank(pheno[,1])-c)/(spl_num-2*c+1))
#pheno[,2]=qnorm((rank(pheno[,2])-c)/(spl_num-2*c+1))

pheno <- sweep(pheno, 2, colMeans(pheno), "-")

task_num=100
test_num_per_task=ceiling(dim(pair_list)[1]/task_num)

task_id =1

	start_idx=(task_id-1)*test_num_per_task+1
	end_idx=task_id*test_num_per_task
	if(end_idx>dim(pair_list)[1]) end_idx=dim(pair_list)[1]

	out_matrix<-matrix(-1.0,(end_idx-start_idx+1),4)

	out_ptr=1
	pre_idx_1=0
	for(pair_idx in start_idx:end_idx)
	{
		
		print(paste(out_ptr," of ",end_idx-start_idx+1,sep=""))

		idx_1=pair_list[pair_idx,1]
		idx_2=pair_list[pair_idx,2]
		gene_name1=geno_info$gene_name[idx_1]
		gene_name2=geno_info$gene_name[idx_2]
		chr1=geno_info$gene_chr[idx_1]
		chr2=geno_info$gene_chr[idx_2]

		if(idx_1!=pre_idx_1) x_A=as.matrix(geno_info$geno_expn[,(geno_info$gene_stnd[idx_1]+1):geno_info$gene_stnd[idx_1+1]])	

		x_B=as.matrix(geno_info$geno_expn[,(geno_info$gene_stnd[idx_2]+1):geno_info$gene_stnd[idx_2+1]])	

		
		out_matrix[out_ptr,1]<- gene_name1
		out_matrix[out_ptr,2]<- gene_name2
		rlt=try(mvFRG_interaction(pheno,x_A,x_B))
		out_matrix[out_ptr,3]<- as.numeric(rlt[1])
		out_matrix[out_ptr,4]<- as.numeric(rlt[2])
		out_ptr=out_ptr+1
		pre_idx_1=idx_1
	}

	write.csv(out_matrix,paste(task_id,"_out.csv",sep=""))


