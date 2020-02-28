#!/usr/bin/perl -w

use strict;
my $name = $ARGV[0];
open(FWD, ">$name.Watson.sam") || die("Error writing to file\n");
open(REV, ">$name.Crick.sam") || die("Error writing to file\n");
while(my $line = <STDIN>){
	my @f = split /\t/, $line;
	if($f[1] eq 16){
		print REV $line;
	}else{
		print FWD $line;
	}
}
close(FWD);
close(REV);
