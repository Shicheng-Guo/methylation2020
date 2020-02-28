#!/usr/bin/perl -w

use strict;

my %hash;

while(my $line = <STDIN>){
	chop($line);
	my @f = split(/\t/, $line);
	next if ($f[3] == 0);
}
