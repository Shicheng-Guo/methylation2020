#!/usr/bin/perl -w
use strict;

#command: /home/ddiep/softwares/bcl2fastq_v1.8.3-build/perl/bin/configureBclToFastq.pl --sample-sheet sample_sheet.csv --input-dir 150331_M00159_0045_000000000-AEUC9/Data/Intensities/BaseCalls --output-dir /media/2TB_Dinh/150331_MiSeq_demultiplex_Dinh --force --use-bases-mask y51,I8,I8 --fastq-cluster-count 0 --flowcell-id BOGUSID

my $lane = 8;
my $project = "Project";
my $i5_file = $ARGV[0];
my $i7_file = $ARGV[1];
my $id_file = $ARGV[2];
my %i5_values;
my %i7_values;

my %ids;

if($id_file){

	open(INFILE, "$id_file") || die("Error opening $id_file\n");
	while(my $line = <INFILE>){
        	chomp($line);
	        my ($sid, $sname) = split "\t", $line;
        	$ids{$sname} = $sid;
	}
	close(INFILE);

	open(IN5, "$i5_file") || die("Error opening $i5_file\n");
	open(IN7, "$i7_file") || die("Error opening $i7_file\n");
	print "FCID,Lane,SampleID,Sample_Ref,Index,Description,Control,Recipe,Operator,SampleProject\n";
	while(my $line = <IN5>){
		chomp($line);
		my ($id, $seq) = split "\t", $line;
		$i5_values{$id} = $seq;
	} 
	while(my $line = <IN7>){
		chomp($line);
		my ($id, $seq) = split "\t", $line;
		$i7_values{$id} = $seq;
	}
	foreach my $i7 (keys %i7_values){
		foreach my $i5 (keys %i5_values){
			my $sname = "$i7-$i5";
			if($ids{$sname}){
				my $sid = $ids{$sname};
				$sname = $sid;
			}
			my ($seq7, $seq5) = ($i7_values{$i7}, $i5_values{$i5});
			print "BOGUSID,$lane,$sname,,$seq7$seq5,,,,,$project\n";
		}
	}
}
else{
	my $in_file = $i5_file;
	open(IN, "$in_file") || die("Error opening $in_file\n");
	print "FCID,Lane,SampleID,Sample_Ref,Index,Description,Control,Recipe,Operator,SampleProject\n";
	while(my $line = <IN>){
		chomp($line);
		my ($id, $seq) = split "\t", $line;
		print "BOGUSID,$lane,$id,,$seq,,,,,$project\n";
	}
}
