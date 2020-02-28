#for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrY chrM chrX
for c in chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrY chrM chrX
do
	sed "s/chr1/$c/g" list_mf > list_mf_${c}
	/home/dinh/scripts/allMethyl2Matrix_DD_06102015.pl list_mf_${c} 10 10 NA freq
done
