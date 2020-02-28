#!/usr/bin/perl -w
use strict;
my $start = $ARGV[0] ? $ARGV[0] : 0;
my $len = $ARGV[1] ? $ARGV[1] : 50;

while(my $line1 = <STDIN>){
	my $line2 = <STDIN>;
	my $line3 = <STDIN>;
	my $line4 = <STDIN>;
	next if (!$line2 || !$line3 || !$line4);
	chop($line1);chop($line2);chop($line3);chop($line4);
	my $actual_len = length($line2)>=($start+$len) ? $len : length($line2)-$start;
	my $seq1 = substr($line2,$start,$actual_len);
	my $seqQual1 = substr($line4,$start,$actual_len);
	print $line1,"\n", $seq1, "\n", $line3, "\n", $seqQual1, "\n";
}
