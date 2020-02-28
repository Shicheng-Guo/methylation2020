email="diep.hue.dinh@gmail.com"
#===Change the following paths===#
cur_dir=`pwd`
scripts_dir="/home/ddiep/scripts/MethylationPipeline/scripts"

# reference files
ref_dir="/oasis/tscc/scratch/ddiep"
ref_fa="$ref_dir/bisHg19/hg19_lambda.fa"
ref_fai="$ref_dir/bisHg19/hg19_lambda.fa.fai"
cpg_list="$ref_dir/bisHg19/hg19_lambda.cpg.positions.txt"

# softwares
samtools="/home/ddiep/softwares/samtools-0.1.18/samtools"
samtools_snp="/home/ddiep/softwares/samtools-0.1.8/samtools"

#===Change the following chromosome names (if other than human)===#
for f in *merged.bam
do
	n=`echo $f | sed 's/.merged.bam//g'`	
	echo $n
	echo "#!/bin/csh" > $n.job
	echo "#PBS -l nodes=1:ppn=2" >> $n.job
	echo "#PBS -l walltime=12:00:00" >> $n.job
	echo "#PBS -o $n.log" >> $n.job
	echo "#PBS -e $n.err" >> $n.job
	echo "#PBS -M $email" >> $n.job
	echo "#PBS -m abe" >> $n.job
	echo "cd $cur_dir" >> $n.job
	echo "$samtools mpileup -B -f $ref_fa $f | $scripts_dir/extractMethyl.pl $cpg_list > $n.methylFreq" >> $n.job
	qsub -q glean $n.job
done
