#!/usr/bin/perl -w
use strict;
my $depletedBase = $ARGV[0];

while(my $line1 = <STDIN>){
	my $line2 = <STDIN>;
	my $line3 = <STDIN>;
	my $line4 = <STDIN>;
	next if (!$line2 || !$line3 || !$line4);
	chomp($line1);chomp($line2);chomp($line3);chomp($line4);
	my $count = ($line2 =~ s/$depletedBase/$depletedBase/g);
	#print $count, "\n";
	next if($count/length($line2) > 0.12);
	print $line1,"\n", $line2, "\n", $line3, "\n", $line4, "\n";
}
