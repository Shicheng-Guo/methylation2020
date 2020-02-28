#!/usr/bin/perl -w
# this script reads in two locus alleles files, calculate the chisqIndep for each comparable locus, and output independent locus records in file A
# combinatorial-entropy-optimization_cgLocusAlleles [file A] [file b]

use strict;
use warnings;
use Statistics::R;
use POSIX;
no warnings 'recursion';

my %allelesTable;
my $lib_size_A = 0;
my $lib_size_B = 0;

printUsage() if(!$ARGV[0]);

my $locus_num_cg;

sub main{
	readCGAllelesTable();
	foreach my $index (keys %allelesTable){
		my @candidates = keys %{$allelesTable{$index}};
		my %alleles;
		my $N_A = 0;
		my $N_B = 0;
		my $foreground_subtract = 0;
		my $background_subtract = 0;
		my (@x, @y);
		foreach my $can (@candidates){
			my ($a, $b) = (0,0);
			$can =~ s/:A//g;
			$can =~ s/:B//g;
			$a = $allelesTable{$index}->{$can.":A"};
			$b = $allelesTable{$index}->{$can.":B"};
			next if(!$a and !$b);
			$a = 0 if(!$a);
			$b = 0 if(!$b);
			$N_A += $a;
			$N_B += $b;
		}
		next if($N_B eq 0 or $N_A eq 0);
		my $normalized_N_A = 200*($N_A/$lib_size_A)/($N_A/$lib_size_A+$N_B/$lib_size_B);
		my $normalized_N_B = 200*($N_B/$lib_size_B)/($N_A/$lib_size_A+$N_B/$lib_size_B);
		#print $normalized_N_A, "\t", $normalized_N_B, "\n";
		foreach my $can (@candidates){
			my ($a, $b) = (0,0);
			$can =~ s/:A//g;
			$can =~ s/:B//g;
			$a = $allelesTable{$index}->{$can.":A"};
			$b = $allelesTable{$index}->{$can.":B"};
			next if(!$a and !$b);
			$a = 0 if(!$a);
			$b = 0 if(!$b);
			$a = $a*$normalized_N_A/$N_A;
			$b = $b*$normalized_N_B/$N_B;
			push(@x, $a);
			push(@y, $b);
			my $N_i_j = $a + $b;
			$foreground_subtract -= logfactorial_gamma($a);
			$foreground_subtract -= logfactorial_gamma($b);
			$background_subtract -= logfactorial_gamma($normalized_N_A*$N_i_j/200);
			$background_subtract -= logfactorial_gamma($normalized_N_B*$N_i_j/200);
		}
		my $output = $foreground_subtract - $background_subtract;
		if(defined($output)){
			print $index, "\t", join(",", @x), "\t", join(",", @y), "\t",$output,"\n"; # if($chi->{"p_value"} < 0.10);		
		}else{
			print $index, "\t", join(",", @x), "\t", join(",", @y), "\tError.Test\n";
		}
	}
	close(IN);
	undef(%allelesTable);
}

sub readCGAllelesTable{
	open(IN, "$ARGV[0]") || die("Error opening $ARGV[0]\n");
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
			next if($a =~ m/[GAN]/);
			$allelesTable{$index}->{$a . ":A"} = $d;
			$lib_size_A += $d;
		}
	}
	close(IN);
	open(IN, "$ARGV[1]") || die("Error opening $ARGV[1]\n");
	while(my $line = <IN>){
		next if($line =~ /indice/);
		chomp($line);
		my @tmp = split /\t/, $line;
		my $index = $tmp[0];
		next if(!$allelesTable{$index});
		my @counts = split ",", $tmp[1];
		my $total_depth = $tmp[2];
		foreach my $value (@counts){
			next if(!$value);
			my ($spacer, $d, $a) = split /(\d+)/, $value;
			$allelesTable{$index}->{$a . ":B"} = $d;
			next if($a =~ m/[GAN]/);
			$lib_size_B += $d;
		}
	}
	close(IN);
}

sub logfactorial{
	my $n = shift;
	return 0 if($n <= 1);
	return log($n)+logfactorial($n-1);
}

sub logfactorial_gamma{
	my $n = shift;
	my $z = $n + 1;
	# only works for positive real numbers.
	my @q = (75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511);
	my $log_factorial = ($z+0.5)*log($z+5.5) - $z - 5.5;
	my $denom_sum = 0;
	for(my $i = 0; $i < 7; $i++){
		$log_factorial -= log($z+$i);
		$denom_sum += $q[$i]*$z^$i;
	}
	$log_factorial += log($denom_sum);
	return $log_factorial;
}

sub printUsage{
	print " Usage: \n";
	print " ./combinatorial-entropy-optimization_cgLocusAlleles.pl [file A] [file B]\n";
	exit 0;
}

main();
