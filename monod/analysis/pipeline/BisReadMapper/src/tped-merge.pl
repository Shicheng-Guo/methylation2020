#!/usr/bin/perl -w

use strict;

my %snp;
my %alleles;
my @names;

while(my $line = <STDIN>){
	chomp($line);
	push(@names, $line);
	open(INFILE, "$line")|| die("ERR opening $line\n");
	while(my $fline = <INFILE>){
		chomp($fline);
		my @f = split /[\ \t]/, $fline;
		if(!$f[4]){
			print $fline,"\n";
			exit 0;
		}
		my $index = join(" ", @f[0...3]);
		$snp{$index}->{$line} = $f[4]. " ". $f[5];
		$alleles{$index}->{$f[4]} = 1;
		$alleles{$index}->{$f[5]} = 1;
	}
	close(INFILE);
}

open(FAM, ">$ARGV[0].tfam") || die("error writing\n");
for(my $i = 0; $i < scalar(@names) ; $i++){
	print FAM "$i ", $names[$i], " 0 0 0 -9\n";
}
close(FAM);

my %diffTable;
my %sumTable;
open(PED, ">$ARGV[0].tped") || die("error writing\n");
foreach my $value(keys %snp){
	my @array = keys %{$alleles{$value}};
	next if(scalar(@array) > 2);
	print PED $value;
	for(my $i = 0; $i < scalar(@names) ; $i++){
		if($snp{$value}->{$names[$i]}){
			print PED " ", $snp{$value}->{$names[$i]};
		}else{ 
			print PED " 0 0";
		}
		for(my $j= 0; $j < scalar(@names); $j++){
			next if(!$snp{$value}->{$names[$i]});
			next if(!$snp{$value}->{$names[$j]});
			$diffTable{$names[$i]}->{$names[$j]}++ if($snp{$value}->{$names[$i]} ne $snp{$value}->{$names[$j]});
			$sumTable{$names[$i]}->{$names[$j]}++;
		}
	}
	print PED "\n";
}
close(PED);

foreach my $value1 (@names){
	print $value1;
	foreach my $value2 (@names){
		if(!$diffTable{$value1}->{$value2}){
			print "\t0";
		}else{
			print "\t", $diffTable{$value1}->{$value2}/$sumTable{$value1}->{$value2};
		}
	}
	print "\n";
}
