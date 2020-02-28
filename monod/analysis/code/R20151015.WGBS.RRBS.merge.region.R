
# test process
scp /home/k4zhang/bin/mergedBam2hapInfo_WGBS_25Jun15.pl shg047@genome-miner.ucsd.edu:/home/shg047/bam/bam2haplo/bam/
scp /home/k4zhang/human_g1k_v37/HsGenome19.CpG.positions.txt shg047@genome-miner.ucsd.edu:/home/shg047/db/hg19


# code 
vim /home/kunzhang/CpgMIP/MONOD/Data/get_methHap_load_matrix_01Oct2015.pl


# bam files


# obtain MHL file for N37
setwd("/home/shg047/bam/bam2haplo")
wgbs<-read.table("/home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined/WGBS_methHap_load_matrix_Oct2015.txt",head=T,sep="\t",row.names=1)
colnames(wgbs)

newdata<-wgbs[,grep("WB",colnames(wgbs))]  # ten columns
write.table(newdata,"WB.mhl.txt",sep="\t",col.names=NA,row.names=T,quote=F)
# N37 bam files
cd /home/k4zhang/my_oasis_tscc/MONOD/N37_WGBS/BAMfiles
cd /media/LTS_33T/WGBS_LTS33/Hg19/Noi_N37_WGBS/BAMfiles
cd /media/LTS_33T/WGBS_LTS33/Hg19/Heyn2013/Re-map/BAMfiles/

# target files
/home/shg047/bam/bam2haplo/target.bed


/home/shg047/bam/bam2haplo/bam



mergedBam2hapInfo.pl target_list_file merged_bam

system()


/home/k4zhang/my_oasis_tscc/MONOD/batch_bam2hapInfoRRBS.pl




cd /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined 


WB_centenarian.all_chrs WB_middle.age.all_chrs  WB_new.born.all_chrs
chr15:89164634-89164825 0.000584795321637427    0.00043010752688172     0.00043010752688172
chr8:11566233-11566251  NA      0.0025  0.0025

WB:WB_centenarian.all_chrs: chr8:11566233-11566251

grep chr8:11566233 WB_centenarian.all_chrs.hapInfo.txt
grep chr8:11566233 WB_middle-age.all_chrs.hapInfo.txt
grep chr8:11566233 WB_new-born.all_chrs.hapInfo.txt

cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined/WB_centenarian.all_chrs.hapInfo.txt ./ 
cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined/WB_middle-age.all_chrs.hapInfo.txt ./ 
cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined/WB_new-born.all_chrs.hapInfo.txt ./ 
  
perl get_methHap_load_matrix_01Oct2015.pl > wb.methhap
grep chr8:11566233 wb.methhap


grep chr8:11566233 WB_centenarian.all_chrs.hapInfo.txt > chr8.WB_centenarian.all_chrs.hapInfo.txt
grep chr8:11566233 WB_middle-age.all_chrs.hapInfo.txt > chr8.WB_middle-age.all_chrs.hapInfo.txt
grep chr8:11566233 WB_new-born.all_chrs.hapInfo.txt > chr8.WB_new-born.all_chrs.hapInfo.txt






# tview check showed centenarian should not be NA
cd /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles
samtools tview  -p chr8:11566233  centenarian.chr8.rmdup.bam /home/shg047/db/hg19.fa
samtools tview  -p chr8:11566233  middle-age.chr8.rmdup.bam /home/shg047/db/hg19.fa
samtools tview  -p chr8:11566233  new-born.chr8.rmdup.bam /home/shg047/db/hg19.fa


cd /media/LTS_33T/WGBS_LTS33/Hg19/Heyn2013/Re-map/BAMfiles
samtools tview  -p chr18:3659912 centenarian.chr18.rmdup.bam /home/shg047/db/hg19/hg19.fa
samtools tview  -p chr18:3659912 middle-age.chr18.rmdup.bam /home/shg047/db/hg19/hg19.fa

# check haplo chr18:3659912-3659992
cd /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles
samtools tview  -p chr18:3659912 centenarian.chr18.rmdup.bam /home/shg047/db/hg19.fa
samtools tview  -p chr18:3659912 middle-age.chr18.rmdup.bam /home/shg047/db/hg19.fa


samtools tview  -p chr10:100206732 centenarian.chr10.rmdup.bam /home/shg047/db/hg19.fa

cp /media/LTS_33T/WGBS_LTS33/Hg19/Heyn2013/Re-map/BAMfiles/centenarian.chr8.rmdup.bam* /home/shg047/bam/bam2haplo/bam
samtools tview  -p chr8:11566233 centenarian.chr8.rmdup.bam /home/shg047/db/hg19.fa
samtools tview  -p chr8:11566233 centenarian.chr8.rmdup.bam /home/shg047/db/hg19/hg19.fa


perl mergedBam2hapInfo_WGBS_25Jun15.pl targetlist.bed centenarian.chr18.rmdup.bam | less -S > resutl2.txt

chr8:11566233-11566251 


cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined/WGBS_methHap_load_matrix_16Oct2015.txt  /home/shg047/monod/nov
cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined/WGBS_methHap_load_matrix_Oct2015.txt /home/shg047/monod/nov


