# postive 
2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.bam
4_On_target_Heart_Casy_N701S504_S4_L001_R1_001.fastq.GFP_sorted.bam

system("samtools view $bam chr6_Ai14GFP:113030794-113030943> $sam.pad");
113029299-113029448
113030794-113030943

perl readStat.pl 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.bam 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.H40.fastq
grep chr 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.H40.fastq.sam | wc -l
grep GFP 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.H40.fastq.sam | wc -l

 
samtools tview 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.bam
chr6_Ai14GFP:113029406

chr6_Ai14GFP:113028721-113028871

vim readStat.pl

bwa index -a bwtsw knockin.fa 
bwa mem  -O 0 knockin.fa 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.H40.fastq > 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.GFP_sorted.H40.fastq.sam

grep 

# background
1_On_target_Liver_Casn_N701S501_S1_L001_R1_001.fastq.Ai14_sorted.bam
3_On_target_Heart_Casn_N701S503_S3_L001_R1_001.fastq.Ai14_sorted.bam

# mouse refrer 113028721-113028871

CTATACTTTCTAGAGAATAGGAACTTCTTAGGGCGGCCGCGGTAT
AAGAAGAAGGCATGAACATGGTTAGCAGAGGCTCTAGA  


Method 2:

system("samtools view $bam chr6_Ai14:113029411-113029433> $sam.pad");

[shg047@tscc-login2 fastq]$ perl readStat.pl 1_On_target_Liver_Casn_N701S501_S1_L001_R1_001.fastq.Ai14_sorted.bam
Reads with >100M: 383189
Reads with Indel: 22744
Reads with SDClip: 5200
Reads of Total: 391453
grep 151M 1_On_target_Liver_Casn_N701S501_S1_L001_R1_001.fastq.Ai14_sorted.pad
151M Number: 355313


perl readStat.pl 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.Ai14_sorted.bam
Reads with >100M: 564402
Reads with Indel: 36884
Reads with SDClip: 83911
Reads of Total: 961929
151M Number: 447894
grep 151M 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.Ai14_sorted.pad


perl readStat.pl 1_On_target_Liver_Casn_N701S501_S1_L001_R1_001.fastq.Ai14_sorted.bam 2_On_target_Liver_Casn.Ai14_H40.fastq


perl readStat.pl 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.Ai14_sorted.bam 2_On_target_Liver_Casy.Ai14_H40.fastq


perl readStat.pl 1_On_target_Liver_Casn_N701S501_S1_L001_R1_001.fastq.Ai14_sorted.bam 1_On_target_Liver_Casn.Ai14_H40.fastq
Reads with >100M: 383189
Reads with Indel: 22744
Reads with SDClip: 5200
Reads of Total: 391453
perl readStat.pl 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq.Ai14_sorted.bam 2_On_target_Liver_Casy.Ai14_H40.fastq
Reads with >100M: 564402
Reads with Indel: 36884
Reads with SDClip: 83911
Reads of Total: 961929

bwa mem  -O 0 knockin.fa 1_On_target_Liver_Casn.Ai14_H40.fastq > 1_On_target_Liver_Casn.Ai14_H40.fastq.sam
bwa mem  -O 0 knockin.fa 2_On_target_Liver_Casy.Ai14_H40.fastq > 2_On_target_Liver_Casy.Ai14_H40.fastq.sam


less -S 2_On_target_Liver_Casy.Ai14_H40.fastq.sam

pcalcute.pl

samtools view -Sb -o 1_On_target_Liver_Casn.Ai14_H40.fastq.bam 1_On_target_Liver_Casn.Ai14_H40.fastq.sam 
samtools view -Sb -o 2_On_target_Liver_Casy.Ai14_H40.fastq.bam 2_On_target_Liver_Casy.Ai14_H40.fastq.sam



perl readStat-2.pl 1_On_target_Liver_Casn.Ai14_H40.fastq.bam
perl readStat-2.pl 2_On_target_Liver_Casy.Ai14_H40.fastq.bam

less -S 1_On_target_Liver_Casn.Ai14_H40.fastq.sam
less -S 2_On_target_Liver_Casy.Ai14_H40.fastq.sam

[shg047@tscc-login2 fastq]$ perl readStat-2.pl 1_On_target_Liver_Casn.Ai14_H40.fastq.bam
Reads with >1M: 351375
Reads with Indel: 24233
Reads with SDClip: 1
Reads of Total: 383189
Indel Ratio:0.0689662041977944
[shg047@tscc-login2 fastq]$ perl readStat-2.pl 2_On_target_Liver_Casy.Ai14_H40.fastq.bam
Reads with >1M: 472225
Reads with Indel: 15629
Reads with SDClip: 1175
Reads of Total: 564402
Indel Ratio:0.033096511196993
[shg047@tscc-login2 fastq]$ samtools view 1_On_target_Liver_Casn_N701S501_S1_L001_R1_001.fastq.Ai14_sorted.bam | less -S


cd /oasis/tscc/scratch/zhl002/kei_050216/fastq
perl readStat.pl 3_On_target_Heart_Casn_N701S503_S3_L001_R1_001.fastq.Ai14_sorted.bam 3_On_target_Heart_Casn.Ai14_H40.fastq 
perl readStat.pl 4_On_target_Heart_Casy_N701S504_S4_L001_R1_001.fastq.Ai14_sorted.bam 4_On_target_Heart_Casy.Ai14_H40.fastq 
bwa mem  -O 0 knockin.fa 1_On_target_Liver_Casn.Ai14_H40.fastq > 3_On_target_Heart_Casn.Ai14_H40.sam 
bwa mem  -O 0 knockin.fa 2_On_target_Liver_Casy.Ai14_H40.fastq > 4_On_target_Heart_Casy.Ai14_H40.sam 
samtools view -Sb -o 3_On_target_Heart_Casn.Ai14_H40.bam 3_On_target_Heart_Casn.Ai14_H40.sam
samtools view -Sb -o 4_On_target_Heart_Casy.Ai14_H40.bam 4_On_target_Heart_Casy.Ai14_H40.sam
perl readStat-2.pl 3_On_target_Heart_Casn.Ai14_H40.bam
perl readStat-2.pl 4_On_target_Heart_Casy.Ai14_H40.bam

cd /oasis/tscc/scratch/zhl002/kei_012616
perl /oasis/tscc/scratch/zhl002/kei_050216/fastq/readStat.pl 1-CAGOn-Cas_S1_L001_R1_001.fastq.Ai14_sorted.bam 1-CAGOn-Cas_S1.Ai14_H40.fastq 
perl /oasis/tscc/scratch/zhl002/kei_050216/fastq/readStat.pl 2-CAGOn-Cas_S2_L001_R1_001.fastq.Ai14_sorted.bam 2-CAGOn-Cas_S1.Ai14_H40.fastq 
bwa mem  -O 0 /oasis/tscc/scratch/zhl002/kei_050216/fastq/knockin.fa 1-CAGOn-Cas_S1.Ai14_H40.fastq > 1-CAGOn-Cas_S1.Ai14_H40.sam 
bwa mem  -O 0 /oasis/tscc/scratch/zhl002/kei_050216/fastq/knockin.fa 2-CAGOn-Cas_S1.Ai14_H40.fastq > 2-CAGOn-Cas_S1.Ai14_H40.sam 
samtools view -Sb -o 1-CAGOn-Cas_S1.Ai14_H40.bam 1-CAGOn-Cas_S1.Ai14_H40.sam
samtools view -Sb -o 2-CAGOn-Cas_S1.Ai14_H40.bam 2-CAGOn-Cas_S1.Ai14_H40.sam
perl readStat-2.pl 1-CAGOn-Cas_S1.Ai14_H40.bam
perl readStat-2.pl 2-CAGOn-Cas_S1.Ai14_H40.bam

less -S 1-CAGOn-Cas_S1.Ai14_H40.sam

cd /oasis/tscc/scratch/zhl002/kei_050216/fastq
perl readStat-2.pl 3_On_target_Heart_Casn.Ai14_H40.bam
perl readStat-2.pl 4_On_target_Heart_Casy.Ai14_H40.bam

cd /oasis/tscc/scratch/zhl002/kei_012616
perl /oasis/tscc/scratch/zhl002/kei_050216/fastq/readStat-2.pl 1-CAGOn-Cas_S1.Ai14_H40.bam
perl /oasis/tscc/scratch/zhl002/kei_050216/fastq/readStat-2.pl 2-CAGOn-Cas_S1.Ai14_H40.bam


cd /oasis/tscc/scratch/zhl002/kei_050216/fastq
grep GAGCAAGGGCGAGGAGC 1_On_target_Liver_Casn_N701S501_S1_L001_R1_001.fastq
grep GAGCAAGGGCGAGGAGC 2_On_target_Liver_Casy_N701S502_S2_L001_R1_001.fastq

cd /oasis/tscc/scratch/zhl002/kei_050216/fastq
grep GAGCAAGGGCGAGGAGC 3_On_target_Heart_Casn.Ai14_H40.fastq 
grep GAGCAAGGGCGAGGAGC 4_On_target_Heart_Casy.Ai14_H40.fastq 

cd /oasis/tscc/scratch/zhl002/kei_012616
grep GAGCAAGGGCGAGGAGC 1-CAGOn-Cas_S1.Ai14_H40.fastq 
grep GAGCAAGGGCGAGGAGC 2-CAGOn-Cas_S1.Ai14_H40.fastq 


chr13:9534035-9534205
chr14:85157462-85157726
chr2:74925547-74925743
chr15:69159421-69159627

perl readStat.pl 3-CAG-OTS1-Cas_S3_L001_R1_001.fastq.Ai14_sorted.bam bedregion.txt 3-CAG-OTS1-Cas_S3.Ai14_sorted
perl readStat.pl 4-CAG-OTS1-Cas_S4_L001_R1_001.fastq.Ai14_sorted.bam bedregion.txt 4-CAG-OTS1-Cas_S4.Ai14_sorted
 OTS1-Cas_S3
 ls *OTS1-Cas_S4*

samtools view 12-CAG-OTS5Cas_S12_L001_R1_001.fastq.Ai14_sorted.bam | less -S 
samtools view 11-CAG-OTS5Cas_S11_L001_R1_001.fastq.Ai14_sorted.bam | less -S 


##‘





