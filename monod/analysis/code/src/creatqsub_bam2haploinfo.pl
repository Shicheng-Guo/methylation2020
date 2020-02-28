#!/usr/bin/perl -w
use strict;
my $bam2hapinfo = "/home/k4zhang/bin/mergedBam2hapInfo.pl";
my $sample_info_file = $ARGV[0];
my %sample_group;

open(INFILE, "$sample_info_file")||die("Error in opening file $sample_info_file.\n");
while(my $line = <INFILE>){
                chop($line);
                my ($id,$bam_file, $target_bed) = split(/\t/, $line);
                next if(!$id || !$bam_file);
                $sample_group{$id}->{"bam_file"}=$bam_file;
                $sample_group{$id}->{"target_bed"}=$target_bed;
}
close(INFILE);

foreach my $id (keys(%sample_group)){
    my $job_file_name = $id . ".job";
    my $status_file = $id.".status";
        my $hapInfo_file = $id.".hapInfo.txt";
        my $curr_dir = `pwd`;
        open(JOB_FILE, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
        print JOB_FILE "#!/bin/csh\n";
        print JOB_FILE "#PBS -q hotel\n";
        print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
        print JOB_FILE "#PBS -l walltime=8:00:00\n";
        print JOB_FILE "#PBS -o ".$id.".log\n";
        print JOB_FILE "#PBS -e ".$id.".err\n";
        print JOB_FILE "#PBS -V\n";
        print JOB_FILE "#PBS -M kun.zhang.ucsd\@gmail.com \n";
        print JOB_FILE "#PBS -m abe\n";
        print JOB_FILE "#PBS -A k4zhang-group\n";
        print JOB_FILE "cd $curr_dir\n";
        my $cmd = "$bam2hapinfo ".$sample_group{$id}->{"target_bed"}." ".$sample_group{$id}->{"bam_file"}." > $hapInfo_file";
        print JOB_FILE "$cmd\n";
        close(JOB_FILE);
        print "Job file $job_file_name created.\n";
        system("qsub $job_file_name");
}
