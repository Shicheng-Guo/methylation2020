#!/usr/bin/perl -w
# this script reads in two locus alleles files, calculate the chisqIndep for each comparable locus, and output independent locus records in file A
# chiSqIndep_cgLocusAlleles [file A] [file b]

use strict;
use warnings;
use Statistics::ChisqIndep;
use POSIX;

my %allelesTable;
my @alleles;

printUsage() if(!$ARGV[0]);
my $out_file = $ARGV[0].".Versus.".$ARGV[1].".Chsq";

my $locus_num_cg;

sub main{
	readCGAllelesTable();
	open(OUT, ">$out_file") || die("Error writing to file $out_file\n");
	print OUT "indices:\tAlleleCounts\tTotalReads\tChisqProbs\n";
	open(IN, "$ARGV[0]") || die("Error reading $ARGV[0]\n");
	while(my $line = <IN>){
		chomp($line);
		my @tmp = split /\t/, $line;
		my $index = $tmp[0];
		next if(!$allelesTable{$index}->{"num"});
		my @counts = split ",", $tmp[1];
		my $total_depth = $tmp[2];
		foreach my $value(@counts){
			next if(!$value || $value !~ /[UM]/);
			my ($spacer, $d, $a) = split /(\d+)/, $value;
			if(!$locus_num_cg){
				$locus_num_cg = length($a);
				generateAlleles();
			}
			$allelesTable{$index}->{$a . ":A"} = $d;
		}
		my @AoA;
		my $count = 0;
		for(my $i = 0; $i < scalar(@alleles); $i++){
			my $a = $allelesTable{$index}->{$alleles[$i].":A"};
			my $b = $allelesTable{$index}->{$alleles[$i].":B"};
			next if(!$a and !$b);
			$a = 0 if(!$a);
			$b = 0 if(!$b);
			push @{ $AoA[0] }, $a;
			push @{ $AoA[1] }, $b;
			$count++;
		}
		next if($count == 1);
		my $chi = new Statistics::ChisqIndep;
		$chi->load_data(\@AoA);
		if($chi->{valid}){
			#print $line, "\n";
			#print $chi->print_summary(), "\n";
			print OUT $line, "\t", sprintf("%8.6f", $chi->{"p_value"}), "\n"; # if($chi->{"p_value"} < 0.10);		
		}else{
			print OUT $line, "\tErrorChiSq\n";
		}
	}
	close(OUT);
	close(IN);
	undef(%allelesTable);
	undef(@alleles);
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
			next if(!$value || $value !~ /[UM]/);
			my ($spacer, $d, $a) = split /(\d+)/, $value;
			$allelesTable{$index}->{$a . ":B"} = $d;
		}
	}
	close(IN);
}

sub generateAlleles{
        @alleles = ("U", "M");
        my $i = 1;
        while($i < $locus_num_cg){
                my @new_alleles;
                foreach my $value (@alleles){
                        push(@new_alleles, $value ."U");
                        push(@new_alleles, $value ."M");
                }
                undef @alleles;
                @alleles = @new_alleles;
                undef @new_alleles;
                $i++;
        }
}

sub printUsage{
	print " Usage: \n";
	print " ./chiSq_cgLocusAlleles.pl [file A] [file B]\n";
	exit 0;
}

main();
