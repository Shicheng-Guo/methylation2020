#===Change the following paths===#
cur_dir=`pwd`
scripts_dir="/home/ddiep/scripts/BisReadMapper/src"
reads_dir="$cur_dir/FASTQ_2"

# reference files
ref_dir="/oasis/tscc/scratch/ddiep"
ref_fai="$ref_dir/bisHg19/hg19_lambda.fa.fai"
template_fwd="$ref_dir/bisHg19/hg19_lambda.bis.CT"
template_rev="$ref_dir/bisHg19/hg19_lambda.bis.GA"

# softwares
bwa="/home/ddiep/softwares/bwa-0.7.5a/bwa"

#=== List FASTQ to process ===#
cd $reads_dir
FILES=`ls *cope.trimmed.fq`
#FILES="s_1_1_ILMN_Indx01.cope.trimmed.fq"
cd $cur_dir

#===Begin===#
for f in ${FILES}
do
	n=`echo $f | sed 's/.fq//g'`
	mkdir $cur_dir/$n
	cd $cur_dir/$n
	#1) Run mapper:
        echo "#!/bin/csh" > $n.job
        echo "#PBS -l nodes=1:ppn=4" >> $n.job
        echo "#PBS -l walltime=14:00:00" >> $n.job
        echo "#PBS -o $n.log" >> $n.job
        echo "#PBS -e $n.err" >> $n.job
        echo "#PBS -V" >> $n.job
        echo "#PBS -M diep.hue.dinh@gmail.com" >> $n.job
        echo "#PBS -m abe" >> $n.job
        echo "#PBS -A k4zhang-group" >> $n.job
        echo "cd /state/partition1/\$USER/\$PBS_JOBID" >> $n.job
        echo "$scripts_dir/BisReadMapper.pl -r $reads_dir/${n}.fq -W $template_fwd -C $template_rev -g $ref_fai -a $bwa -b 64 -p 4 -n $n > $n.status" >> $n.job
	echo "cp *sam $cur_dir/$n" >> $n.job
	echo "cp *status $cur_dir/$n" >> $n.job
	#echo "cd $cur_dir/$n" >> $n.job
	echo "ls *sorted.sam > list_sams" >> $n.job
        echo "$scripts_dir/BamExtractor.pl -i list_sams -s $cur_dir/list_paths_Hg19 -o $n -r none -v no -b no -p no -d 5" >> $n.job
	echo "cp *bam $cur_dir/$n" >> $n.job
        echo "cp *methylFreq $cur_dir/$n" >> $n.job
	#qsub -q hotel $n.job
done
#===End===#
