#!/usr/bin/perl -w

my %hash;

while(my $line = <STDIN>){
	chomp($line);
	$hash{$line}++;
}

foreach my $value(keys %hash){
	print $value, "\t", $hash{$value}, "\n";
}
