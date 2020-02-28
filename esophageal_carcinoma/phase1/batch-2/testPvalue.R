combineTestPvalsMeth <- function(pvalues) {
    # use Fisher's classical method for combining p-values in the case of independence
    tcombined = sum(-2*log(pvalues))
    return(pchisq(tcombined,2*length(pvalues),lower.tail=FALSE))
}

methy_ttest=function(data){
  type01=which(data[1,]=="01")
  type11=which(data[1,]=="11")
  test.result=apply(data[-1,],1,function(x) {a=t.test(as.numeric(x[type01]),as.numeric(x[type11])); b=as.numeric(a$p.value); c=as.numeric(data.frame(a$estimate)[1:2,1]);  e=ifelse(as.numeric(c[1])>as.numeric(c[2]),"Hyper","Hypo"); d=append(append(b,c),e); return(d)  })
  rownames(test.result)=c("pvalue","cancer_avg","normal_avg","HyperorHypo")
  test.result=t(test.result)
  test.result=as.data.frame(test.result)
  rownames(test.result)=rownames(data)[-1]
  test.result[,1:3]=apply(test.result[,1:3],2,function(x) as.numeric(x))
  return(test.result)
  }



