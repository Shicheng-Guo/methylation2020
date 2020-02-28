#!/usr/bin/perl -w

use strict;

my $fraction = $ARGV[0];
my $randomNumber = int($fraction*1000);

#chr22:45103758-45103899 CCCCCC  1       45103765,45103814,45103836,45103846,45103849,45103857
while(my $line = <STDIN>){
	chomp($line);
	my ($probeID, $hapString, $count, $cpgPos) = split /\t/, $line;
	my $i = $count;
	my $keep = 0;
	do{
		my $number = int(rand()*1000);
		$keep++ if($number <= $randomNumber);
		$i--;
	}while($i > 0);
	if($keep > 0){
		print "$probeID\t$hapString\t$keep\t$cpgPos\n";
	}
}

