#!/usr/bin/perl -w
# this script reads in two locus alleles files, calculate the chisqIndep for each comparable locus, and output independent locus records in file A
# kolmogorov-smirnov_cgLocusAlleles [file A] [file b]

use strict;
use warnings;
use Statistics::R;
use POSIX;

my %allelesTable;

printUsage() if(!$ARGV[0]);

my $locus_num_cg;

sub main{
	readCGAllelesTable();
	open(IN, "$ARGV[0]") || die("Error reading $ARGV[0]\n");
	while(my $line = <IN>){
		chomp($line);
		my @tmp = split /\t/, $line;
		my $index = $tmp[0];
		next if(!$allelesTable{$index}->{"num"});
		my @counts = split ",", $tmp[1];
		my $total_depth = $tmp[2];
		foreach my $value(@counts){
			next if(!$value);
			my ($spacer, $d, $a) = split /(\d+)/, $value;
			#print "$spacer,$d,$a\n";
			if(!$locus_num_cg){
				$locus_num_cg = length($a);
			}
			$allelesTable{$index}->{$a . ":A"} = $d;
		}
		my @x;
		my @y;
		my $count = 0;
		my @candidates = keys %{$allelesTable{$index}};
		my %alleles;
		foreach my $can (@candidates){
			$can =~ s/:A//g;
			$can =~ s/:B//g;
			$alleles{$can} = 1;
		}
		foreach my $hap (keys %alleles){
			my $a = $allelesTable{$index}->{$hap.":A"};
			my $b = $allelesTable{$index}->{$hap.":B"};
			next if(!$a and !$b);
			$a = 0 if(!$a);
			$b = 0 if(!$b);
			push @x, $a;
			push @y, $b;
			$count++;
		}
		next if($count < 3);
		my $R = Statistics::R->new();
		$R->startR;
		$R->set('x', \@x );
		$R->set('y', \@y );
		$R->run(q`z <- ks.test(x,y)`);
		my $output = $R->get('z$p.value');
		$R->stopR;
		if($output){
			#print $line, "\n";
			#print $chi->print_summary(), "\n";
			print $index, "\t", join(",", @x), "\t", join(",", @y), "\t",$output,"\n"; # if($chi->{"p_value"} < 0.10);		
		}else{
			print $index, "\t", join(",", @x), "\t", join(",", @y), "\tErrorKS.Test\n";
		}
	}
	close(IN);
	undef(%allelesTable);
}

sub readCGAllelesTable{
	open(IN, "$ARGV[1]") || die("Error opening $ARGV[1]\n");
	while(my $line = <IN>){
		next if($line =~ /indice/);
		chomp($line);
		my @tmp = split /\t/, $line;
		my $index = $tmp[0];	
		my @counts = split ",", $tmp[1];
		my $total_depth = $tmp[2];
		$allelesTable{$index}->{"num"} = 1;
		foreach my $value (@counts){
			next if(!$value);
			my ($spacer, $d, $a) = split /(\d+)/, $value;
			$allelesTable{$index}->{$a . ":B"} = $d;
		}
	}
	close(IN);
}

sub printUsage{
	print " Usage: \n";
	print " ./kolmogorov-smirnov_cgLocusAlleles.pl [file A] [file B]\n";
	exit 0;
}

main();
