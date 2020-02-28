#!/usr/bin/perl -w
use strict;

my $target_bed = $ARGV[0];

my %cntTable;
my %foundTable;
open(INFILE, "$target_bed") || die("Error in opening targetFile!\n");
while(my $line = <INFILE>){
	chop($line);
	next if($line !~ /^chr[0-9XY]/);
	my ($chr,$start,$end) = split(/[\t ]+/, $line);
	$chr =~ s/chr//;
	if($start > $end){
		my $temp = $end;
		$end = $start;
		$start = $end;
	}
	for(my $i = $start; $i<=$end; $i++){
		$foundTable{$chr .":". $i} = 1;
	}
}
close(INFILE);

my $total_target_size = scalar(keys %foundTable);

my @datasets;

my $file = $ARGV[1];	
open(INFILE, "$file") || die("Error in opening $file!\n");
	$cntTable{$file}->{"total"} = 0;
	while(my $line = <INFILE>){
		my @fields = split(/\t/, $line);
		my $chr = $fields[0];
		$chr =~ s/chr//;
		$cntTable{$file}->{"total"} +=$fields[3];
		if($foundTable{$chr. ":".$fields[1]}){
			$cntTable{$file}->{$chr.":".$fields[1]} += $fields[3];
		}
	}
	push(@datasets, $file);
close(INFILE);

foreach my $file(@datasets){
	my $total_target_bases = 0;
	my $total_target_covered = 0;
	foreach my $index(keys %foundTable){
		next if(!$cntTable{$file}->{$index});
		my $count = $cntTable{$file}->{$index};
		$total_target_bases += $count;
		$total_target_covered ++ if($count > 0);
	}	
	my $total_bases = $cntTable{$file}->{"total"};
	print "Size of target = $total_target_size\n";
	print "Total mappable bases = $total_bases\n";
	print "Total target bases = $total_target_bases (", sprintf("%4.3f", $total_target_bases/$total_bases), ") \n";
	print "Total target covered 1X = $total_target_covered (", sprintf("%4.3f", $total_target_covered/$total_target_size), ") \n";
}
