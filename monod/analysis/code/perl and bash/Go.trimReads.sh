reads_dir="/oasis/tscc/scratch/ddiep/Working/140722_RRBS/FASTQ/140627_MISEQ"
cur_dir=`pwd`
cope="/home/ddiep/softwares/cope-src-v1.1.3/src/cope -o connect.fq -2 left1.fq -3 left2.fq -m 0"
#THIS ONE for PE 
trim="/home/ddiep/softwares/trim_galore_latest/trim_galore --phred64 --rrbs --paired --dont_gzip -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

FILES="NP-RRBS-NC-P-1ng-p2-Jun20-A013_S3_L001_R1_001.fastq.gz NP-RRBS-NC-P-1ng-p2-Jun20-A020_S6_L001_R1_001.fastq.gz NP-RRBS-NC-P-1ng-p2-Jun20-A018_S4_L001_R1_001.fastq.gz NP-RRBS-PC-P-1ng-p2-Jun20-A010_S1_L001_R1_001.fastq.gz NP-RRBS-NC-P-1ng-p2-Jun20-A019_S5_L001_R1_001.fastq.gz NP-RRBS-PC-P-1ng-p2-Jun20-A011_S2_L001_R1_001.fastq.gz"

cd $cur_dir

for f in ${FILES}
do
	n=`echo $f | sed 's/_L001_R1_001.fastq.gz//g'`
	g=`echo $f | sed 's/R1/R2/g'`
        echo "#!/bin/csh" > $n.job
        echo "#PBS -l nodes=1:ppn=2" >> $n.job
        echo "#PBS -l walltime=3:00:00" >> $n.job
        echo "#PBS -o $n.log" >> $n.job
        echo "#PBS -e $n.err" >> $n.job
        echo "#PBS -V" >> $n.job
        echo "#PBS -M diep.hue.dinh@gmail.com" >> $n.job
        echo "#PBS -m abe" >> $n.job
        echo "#PBS -A k4zhang-group" >> $n.job
        echo "cd /state/partition1/\$USER/\$PBS_JOBID" >> $n.job	
	echo "$trim $reads_dir/$f $reads_dir/$g" >> $n.job
	v1=`echo $f | sed 's/.fastq.gz//g'`
	v2=`echo $g | sed 's/.gastq.gz//g'`
	echo "$cope -a ${v1}_val_1.fq -b ${v2}_val_2.fq > $n.cope.log.txt" >> $n.job
	echo "cat connect.fq left1.fq left2.fq > $n.cope.trimmed.fq" >> $n.job
	echo "cp *fq $cur_dir/" >> $n.job
	echo "cp *txt $cur_dir/" >> $n.job
	qsub -q hotel $n.job
done


