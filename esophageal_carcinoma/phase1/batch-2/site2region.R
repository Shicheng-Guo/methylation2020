annotation = read.csv("f:/Esphageal Cancer/GPL13534_HumanMethylation450_15017482_v.1.1.csv",header=T,sep="\t")
sites2region = function(dmtp, window = 4){
  chr =c(1:22)
  sites2region_list = list()
  list_names <<-c()
  seq = match(rownames(dmtp), annotation[,1])
  dmtp = data.frame(dmtp, annotation[seq,11:12])
  dmtp$MAPINFO = as.numeric(dmtp$MAPINFO)
  sites2region_list <<- list()
  list_names <<- c()
  sapply(1:length(chr), FUN =function(i){
    sub_dmtp = dmtp[which(dmtp$CHR == i),]
    ordered = order(sub_dmtp$MAPINFO)
    seq = match( rownames(sub_dmtp)[ordered], rownames(dmtp))
    n_cgsite = dim(sub_dmtp)[1]
    N_groups = n_cgsite - window + 1
    name_1 = paste("CHR",i,sep="")
    name_2 = paste("Group",1:N_groups,sep="")
    names = paste(name_1, name_2, sep="_") 
    list_names <<-c(list_names, names)
    sapply(1:N_groups, FUN =function(j){    
        sites2region_list[[length(sites2region_list)+1]] <<- seq[j:(j+window-1)] 
    })
  })
  names(sites2region_list) = list_names
  return (sites2region_list)
}
