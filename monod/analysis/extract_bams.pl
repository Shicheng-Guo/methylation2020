#!/usr/local/bin/perl -w
use strict;

my $bam_file_list = $ARGV[0];
my $target_bed_file = $ARGV[1];
$target_bed_file= "/home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/Primary_tumor_ALL.genomecov.RD50_80UP.merged.bed" if(!$target_bed_file);
my %bamInfoTable;
my $bedtool_exe = "/home/kunzhang/softwares/bedtools-2.17.0/bin/bedtools";
open(INFILE, "$bam_file_list")||die("Error in opening file $bam_file_list\n");
while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		if(-e $fields[1] && $fields[0]){
			$bamInfoTable{$fields[0]}=$fields[1];
		}
}
close(INFILE);

foreach my $id (sort keys(%bamInfoTable)){
	my $extracted_bam = $id . ".bam";
	my $cmd = "$bedtool_exe intersect -abam $bamInfoTable{$id} -b $target_bed_file > $extracted_bam";
	system($cmd);
	$cmd = "/home/kunzhang/bin/samtools index $extracted_bam";
	system($cmd);
}

