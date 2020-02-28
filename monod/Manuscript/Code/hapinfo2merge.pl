#!/usr/bin/perl -w

# A perl script to merge SRR hapinfo files to single one SRX file.
# Contact: Shihcheng.Guo@Gmail.com
# Go to http://sra.dnanexus.com/studies/SRP028600/samples
# Select SRS and Click Related RUNS then get the Table as the input
# Version 1.3
# Date: 2017-05-06 (publish enough)


use strict;
use warnings;
use Cwd;

use strict;
use warnings;
use Cwd;

my $file=shift @ARGV;
my %SRA;
open F,$file;
while(<F>){
chomp;
if(/(SRR\d+)/){
        my $SRR=$1;
        if(/(SRX\d+)/){
                my $SRX=$1;
                print "$SRR\t$SRX\n";
                push @{$SRA{$SRX}},$SRR;
                }
        }
}
close F;

foreach my $SRX(sort keys %SRA){
        open OUT,">$SRX.hapInfo.txt" if scalar(@{$SRA{$SRX}})>=1;
        my %data;
        print "\t$SRX\n";
        foreach my $SRR (@{$SRA{$SRX}}){
        	print "\tt$SRR\n";
        	open F1,"$SRR.hapInfo.txt" || warn "cannot find or connect $SRR.hapInfo.txt\n";
        	while(<F1>){
				chomp;
				my ($mhb,$haptype,$number,$pos)=split /\t/;
		        next if length($haptype)<2;
				my $key="$mhb\_$haptype\_$pos";
				$data{$key}+=$number
			}
			close(F1)
        }

		foreach my $key(sort keys %data){
			my ($mhb,$haptype,$pos)=split /_/,$key;
			print OUT "$mhb\t$haptype\t$data{$key}\t$pos\n"; 
		}
		close OUT;
		print "$SRX Haplotype file merging are completed!\n";
}


