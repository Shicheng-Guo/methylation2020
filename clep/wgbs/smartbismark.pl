#!/usr/bin/perl -w

# bismark to alignment pair-end targeted methylation sequencing from multiplex pcr (targeted bisulfite sequencing)
# Contact: Shicheng Guo
# Version 1.4
# Update: 10/02/2019
# smartbismark.pl 
# USAGE: perl smartbismark.pl --input SraRunTable.txt --genome hg19 --phred=33 --server MCRI --queue shortq --BismarkRefereDb --submit 

use strict;
use Cwd;
use Getopt::Long;
my $smartbismark_version="v0.16.3_dev";
my $dir=getcwd;
my ($input,$genome,$server,$BismarkRefereDb,$queue,$phred,$help,$submit,$nodes,$ppn,$walltime,$multicore)=process_command_line();
my $sed=directory_build();

open SRR,$input || die("Error in opening file $input.\n");   ;
my %SRR; my %SRA;

while(<SRR>){
    chomp;
    next if /^\s+$/;
	next if /sample_accession|BioProject/i;
	my $SRR;
	if(/(SRR\d+)/){
	$SRR=$1;
	if(/(SRX\d+)/){
		my $SRX=$1;
		$SRR{$SRR}=$SRR;
		push @{$SRA{$SRX}},$SRR;
		}
	}elsif(/(.+)_R1/){
		$SRR=$1;
		my $SRX=$SRR;
		push @{$SRA{$SRX}},$SRR;
	}		
print "$SRR\n";
my @read = glob("$SRR\_*fastq*gz");
my $chrLenhg19="~/hpc/db/hg19/hg19.chrom.sizes";
my $cpgPoshg19="~/hpc/db/hg19/HsGenome19.CpG.positions.txt"; 
my $curr_dir = $dir;
if(scalar(@read) eq 2){
my($sample,undef)=split /_1.fastq.gz|_R1.fastq.gz/,$read[0]; 	
my($sample1,undef)=split /.fastq.gz/,$read[0];
my($sample2,undef)=split /.fastq.gz/,$read[1];

my $job_file_name = "$SRR.pbs";
open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");   
    
print OUT "#!/bin/bash\n";
print OUT "#PBS -N $sample\n";
print OUT "#PBS -q $queue\n";
print OUT "#PBS -l nodes=$nodes:ppn=$ppn\n";
print OUT "#PBS -M Guo.Shicheng\@marshfieldresearch.org\n";
print OUT "#PBS -m abe\n";
print OUT "cd $curr_dir\n";    

my $extractor="bismark_methylation_extractor --no_overlap --multicore $multicore --merge_non_CpG --bedGraph --cutoff 10 --ignore 1 --buffer_size 4G";
# my $bismark="bismark --bowtie2 --non_directional --multicore $multicore --fastq -N 1"; 
my $bismark="bismark --bowtie2 --multicore $multicore --fastq -N 1"; 
my $coverage2cytosine="coverage2cytosine --merge_CpG --gzip --genome_folder";
    
print OUT "# fastq-dump --skip-technical --split-files --gzip $SRR\n";
print OUT "# trim_galore --paired --stringency 3 --phred$phred --clip_R1 2 --clip_R2 2 --fastqc --illumina $sample1\.fastq.gz $sample2\.fastq.gz --output_dir ../fastq_trim\n";
print OUT "# $bismark --phred$phred-quals $BismarkRefereDb -1 ../fastq_trim/$sample1\_val_1.fq.gz -2 ../fastq_trim/$sample2\_val_2.fq.gz -o ../bam\n";
print OUT "# filter_non_conversion --paired ../bam/$sample1\_val_1_bismark_bt2_pe.bam\n";
print OUT "# deduplicate_bismark --bam ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam\n";   
print OUT "# $extractor --comprehensive --output /gpfs/home/guosa/PRJNA528657/methyfreq  ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam\n";
print OUT "$coverage2cytosine $BismarkRefereDb -o /gpfs/home/guosa/PRJNA528657/bed/$sample1.mergeCpG.bed /gpfs/home/guosa/PRJNA528657/methyfreq/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz\n";
print OUT "# samtools sort -@ 8 ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam -o ../sortbam/$sample\_bismark_bt2_pe.sortc.bam\n";
print OUT "# samtools index ../sortbam/$sample\_bismark_bt2_pe.sortc.bam\n";
print OUT "# cd ../sortbam\n";
print OUT "# perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/oasis/db/hg19/hg19.cut10k.bed > saminfo.txt\n";
print OUT "# perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt submit bismark $chrLenhg19 $cpgPoshg19\n";

}elsif(scalar(@read) == 1){
	my($sample,undef)=split /_1.fastq.gz|_R1.fastq.gz/,$read[0]; 	
	my($sample1,undef)=split /.fastq.gz/,$read[0];
	my $job_file_name = "$SRR.pbs";
	open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");  
	my $extractor="bismark_methylation_extractor --single-end --merge_non_CpG --bedGraph --cutoff 5 --ignore 1 --buffer_size 4G";
	print OUT "#!/bin/csh\n";
	print OUT "#PBS -N $sample\n";
	print OUT "#PBS -q $queue\n";
	print OUT "#PBS -l nodes=$nodes:ppn=$ppn\n";
	print OUT "#PBS -l walltime=$walltime\n";
	print OUT "#PBS -o $sample.log\n";
	print OUT "#PBS -e $sample.err\n";
	print OUT "#PBS -V\n";
	print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
	print OUT "#PBS -m abe\n";
	print OUT "#PBS -A k4zhang-group\n";
	print OUT "cd $curr_dir\n";  
	print OUT "fastq-dump --skip-technical --split-files --gzip $sample\n";
	print OUT "trim_galore --phred$phred --fastqc --illumina $sample1.fastq.gz --output_dir ../fastq_trim\n";
	print OUT "bismark --bowtie2 --phred$phred-quals --fastq -N 1 --multicore $multicore $BismarkRefereDb ../fastq_trim/$sample1\_trimmed.fq.gz -o ../bam\n";  
	print OUT "filter_non_conversion --single ../bam/$sample1\_trimmed_bismark_bt2.bam\n"; 
	print OUT "deduplicate_bismark --bam ../bam/$sample1\_trimmed_bismark_bt2.nonCG_filtered.bam\n";
	print OUT "$extractor --comprehensive --output ../methyfreq  ../bam/$sample1\_trimmed_bismark_bt2.nonCG_filtered.deduplicated.bam\n";
	print OUT "samtools sort ../bam/$sample1\_trimmed_bismark_bt2.nonCG_filtered.deduplicated.bam -o ../sortbam/$sample.bismark_bt2_se.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample.bismark_bt2_se.sort.bam\n";
	print OUT "perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/oasis/db/hg19/hg19.cut10k.bed > saminfo.txt\n";
	print OUT "perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt submit bismark $chrLenhg19 $cpgPoshg19\n";
}
close(OUT);
 
if($submit eq 'submit'){
	system("qsub $SRR.pbs");
}else{
	print "qsub $SRR.pbs\n";
}
}

sub process_command_line{
	my $help;
	my $version;
	my $input;
	my $bismark_version;
	my $submit;
	my %walltime;
	my %ppn;
	my %multicore;
	my $nodes;
    my $BismarkRefereDb;
	my $phred;
    
	my $command_line=GetOptions ( 
                                  "input=s"    		    => \$input,
                                  "genome=s"   			=> \$genome,
                                  "server=s"   			=> \$server,
	                              "queue=s"   			=> \$queue,                                                                    
                                  "help"      			=> \$help,
	                              "submit=s"   			=> \$submit,
	                              "BismarkRefereDb=s"   => \$BismarkRefereDb,
	                              "phred=s"             => \$phred,
	                             );
	
    unless (defined $genome){
    warn "\n\n\t[Error]: Please respecify essential command line options!\n";
    print_helpfile();
	}
    #################################################################################################
    ##################### SmartBismark Version and Usage (Version and Usage) ########################
    #################################################################################################
    if($help){
    print_helpfile();
    exit;
    }   
    if(! defined $phred){
	 $phred=33;
	}
    #################################################################################################
    ##################### Assemble Bismark Reference (hg19,hg38,mm9,mm10) ##########################
    #################################################################################################
    if(!defined $BismarkRefereDb){
    if($server eq "TSCC"){
    	if($genome eq "hg19"){
    	$BismarkRefereDb="/home/shg047/db/hg19/bismark/";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome eq "hg38"){
		$BismarkRefereDb="/home/shg047/db/hg38/bismark/";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome eq "mm9"){
		$BismarkRefereDb="/home/shg047/oasis/db/mm9/bismark";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome eq "mm10"){
		$BismarkRefereDb="/home/shg047/oasis/db/mm10/bismark";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}else{
		warn("Please assign genome version (in TSCC)to the script: hg19? hg38? mm9? mm10?");	
		}
    }elsif($server eq "GM"){
		if($genome eq "hg19"){
		$BismarkRefereDb="/media/Home_Raid1/shg047/db/hg19/bismark";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome=="hg38"){
		$BismarkRefereDb="/home/shg047/db/hg38/bismark/";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome eq "mm9"){
		$BismarkRefereDb="/home/shg047/db/mm9/bismark/";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome eq "mm10"){
    	        $BismarkRefereDb="/home/shg047/db/mm10/bismark/";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}else{
		warn("Please assign genome version (in Genome-miner)to the script: hg19? hg38? mm9? mm10?");	
		}
	}elsif($server eq "MCRI"){
		if($genome eq "hg19"){
		$BismarkRefereDb="~/hpc/db/hg19/bismark";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome=="hg38"){
		$BismarkRefereDb="~/hpc/db/hg38/bismark/";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome eq "mm9"){
		$BismarkRefereDb="~/hpc/db/mm9/bismark/";
		print "Bismark alignment reference: $BismarkRefereDb\n";
		}elsif($genome eq "mm10"){
    	        $BismarkRefereDb="~/hpc/db/mm10/bismark/";
				print "Bismark alignment reference: $BismarkRefereDb\n";
		}else{
				warn("Please assign genome version (in Genome-miner)to the script: hg19? hg38? mm9? mm10?");	
		}
	}else{
        print "Please check the server name or build corresponding reference for server:$server\n";
	}
	}
    #################################################################################################
    ##################### Assemble PBS Paramters (nodes, ppn, walltime multicore) ###################
    #################################################################################################
    warn "\nYou didn't assign --queue for the script, the default setting: multicore=2 and ppn=6 will be applied!\n\n" if ! defined $queue; 
    $queue="hotel" if ! defined $queue;
    %walltime=(
    hotel   => "168:00:00",
    condo   => "8:00:00",
    pdafm   => "72:00:00",
    glean   => "72:00:00",
    default => "168:00:00",
    shortq => "240:00:00",
    longq => "480:00:00",
    );
    %ppn=(
    hotel   => "16",
    pdafm   => "32",
    glean   => "16",
    condo   => "16",
    default => "6",
    shortq  => "1",
    longq   => "12",
    );
    %multicore=(
    hotel   => "6",
    pdafm   => "12",
    glean   => "6",
    condo   => "6",
    default => "2",
    shortq  => "1",
    longq   => "4",
    );
    warn "Queue: $queue is not found in this server, please check the queue name\n" if ! defined $ppn{$queue}; 
    warn "Queue: $multicore is not found in this server, please check the queue name\n" if ! defined $multicore{$queue}; 
    $nodes=1;
    $ppn=$ppn{$queue};
    $walltime=$walltime{$queue};
    $multicore=$multicore{$queue};
    return($input,$genome,$server,$BismarkRefereDb,$queue,$phred,$help,$submit,$nodes,$ppn,$walltime,$multicore);
}


sub print_helpfile{
	print<<"HOW_TO";

	USAGE: smartbismark --input saminfo.txt --genome hg19 --server MCRI --submit
	
	ARGUMENTS:
	
	Last edited on 15 Apirl 2017 <contact: Shihcheng.Guo\@gmail.com>.
	
	--input    Configure files including SRR and SRX information. Script will extract SRR and SRX
		   information automatically and then download the corresponding SRR files and then do
	           the further bismark alignment, sort and hapinfo calling. Single-end or pair-end fas
		   tq file list. For single-end fastq, one file in each line. For pair-end fastq files,
                   paired fastq files should be listed in one line with TAB. 
	
	--genome   hg19, hg38, mm9, mm10. The option is used to determined which bismark alignment ref-
	           erence will be chosed. if you have new reference rather than huamn and mouse, please
	           creat the methylation alignment reference advanced and you can contact me to do that. 
	
	--server   TSCC, GM(GenomeMiner),MCRI. Combined with genome, --server and --genome will provide the 
	           location of the alignment reference for bismark. 
	           
	--submit   Submit pbs job or not. SmartBismark will creat pbs job files for each fastq file and
	           defaulty taken the system is PBS system. if --submit="submit", then PBS job will be 
	           submitted and PBS ID will be printed in the STANDOUT.    
	--queue    TSCC queue: hotel, glean, pdafm, condo
			   
HOW_TO
exit;
}

sub print_bismark_version{
	
my $bismark_version=bismark_version();	
	
	print<<"VERSION";
          SmartBismark - Smart Bisulfite Mapper and Methylation Caller Assembly Tools.

                       SmartBismark Version: $smartbismark_version
                          Bismark Version: $bismark_version
                           
        Copyright 2010-15 Shicheng Guo, University of Wisconsin at Madison (UW-Madison)

VERSION
    exit;
}


sub directory_build{
chdir getcwd;
print "===================================================================================\n";
mkdir "../fastq_trim" if ! -e "../fastq_trim" || print  "fastq_trim     Building Succeed!  <Trimed Fastq>     stored here\n";
mkdir "../bam" if ! -e "../bam"               || print  "bam            Building Succeed!    <Bam Files>      stored here\n";
mkdir "../bedgraph" if ! -e "../bedgraph"     || print  "bedgraph       Building Succeed!  <Bedgraph Files>   stored here\n";
mkdir "../sortbam" if ! -e "../sortbam"       || print  "sortbam        Building Succeed!  <SortBam Files>    stored here\n";
mkdir "../methyfreq" if ! -e "../methyfreq"   || print  "methyfreq      Building Succeed!  <MethyFreq Files>  stored here\n";
mkdir "../bedgraph" if ! -e "../bedgraph"     || print  "bedgraph       Building Succeed!  <BedGraph Files>   stored here\n";
mkdir "../bw" if ! -e "../bw"                 || print  "bigwig         Building Succeed!  <BigWig Files>     stored here\n";
mkdir "../hapinfo/" if ! -e "../hapinfo/"     || print  "hapinfo        Building Succeed!  <hapinfo Files>    stored here\n";
mkdir "../bed/" if ! -e "../bed/"             || print  "bed            Building Succeed!  <bed Files>        stored here\n";
print "===================================================================================\n";
}
sub bismark_version{
	my @version=`bismark --version`;
	my $bismark_version;
	foreach my $line(@version){
		if($line=~/version/i){
			$bismark_version=$line;
		}
	}
	return($bismark_version);
}
