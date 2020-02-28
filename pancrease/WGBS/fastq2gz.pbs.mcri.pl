#!/usr/bin/perl -w
use strict;
use Cwd;
my $dir=getcwd;
my $submit=shift @ARGV;
my @fastq=glob("*fq *fastq");
foreach my $fastq(@fastq){
        my $job_file_name = $fastq.".job";
        open JOB_FILE, ">$job_file_name" || die("Error in opening file $job_file_name.\n");
        print JOB_FILE "#!/bin/csh\n";
        print JOB_FILE "#PBS -q shortq\n";
        print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
        print JOB_FILE "#PBS -o ".$fastq.".log\n";
        print JOB_FILE "#PBS -e ".$fastq.".err\n";
        print JOB_FILE "#PBS -V\n";
        print JOB_FILE "#PBS -M Guo.Shicheng\@marshfieldresearch.org\n";
        print JOB_FILE "#PBS -m abe\n";
        print JOB_FILE "cd $dir\n";
        print JOB_FILE "gzip $fastq\n";
        close(JOB_FILE);
        print "$job_file_name created\n";
        if($submit eq "submit"){
        system("qsub $job_file_name");
        }
}
