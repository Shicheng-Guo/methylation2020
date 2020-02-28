#!/usr/bin/perl -w

use strict;
use List::Util qw(shuffle);

my $total_num_files = $ARGV[0];
my $out_prefix = $ARGV[1];
my @numbers;
my %fileHandles;

for(my $i = 1; $i <= $total_num_files; $i++){
	my $n = sprintf("%03d", $i);
	open($fileHandles{$n}, ">$out_prefix.$n.haploInfo.txt");
	push(@numbers, $n);
}

#chr22:45103758-45103899 CCCCCC  1       45103765,45103814,45103836,45103846,45103849,45103857
while(my $line = <STDIN>){
	chomp($line);
	my ($probeID, $hapString, $count, $cpgPos) = split /\t/, $line;
	my $i = $count;
	my %keep;
	do{
		my @shuffledNumbers = shuffle(@numbers);
		my $n = pop(@shuffledNumbers);
		$keep{$n}++;
		$i--;
	}while($i > 0);
	foreach my $n (keys %keep){
		print {$fileHandles{$n}} "$probeID\t$hapString\t", $keep{$n}, "\t$cpgPos\n";
	}
}

foreach my $n (keys %fileHandles){
	close($fileHandles{$n});
}

