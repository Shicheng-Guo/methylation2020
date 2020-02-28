#!/usr/bin/env perl
#
#    Copyright (C) 2014, 2015 Institute for Genomic Medicine (IGM).
#    Portions copyright (C) 2015 Unversity of California, San Diego.
#
#    Author: Kun Zhang <k4zhang@ucsd.edu> and Shicheng Guo <scguo@ucsd.edu>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

#!/usr/bin/perl -w
use strict;
use Cwd;
# my $input_dir=shift @ARGV;
# my $output_dir=shift @ARGV;
my $input_dir="/home/k4zhang/my_oasis_tscc/MONOD/N37_WGBS/BAMfiles";
my $output_dir="/oasis/tscc/scratch/shg047/N37";
# die "USAGE: mergeBam.pl input_dir output_dir \n" if scalar @ARGV <1;

chdir $input_dir || die "can't enter to $input_dir\n";

print "cd /home/k4zhang/my_oasis_tscc/MONOD/N37_WGBS/BAMfiles\n";
my %sample;
my @bam=glob("*.bam");
foreach my $bam(@bam){
my ($sam,$chr,undef)=split /\./,$bam;
push (@{$sample{$sam}},$bam);
}

mkdir $output_dir if (!-e $output_dir);
foreach my $sam(sort keys %sample){
	print "samtools view -H $sample{$sam}[0] > $output_dir/$sam.bam.header\n";
	print "samtools merge -h $output_dir/$sam.bam.header $output_dir/$sam.bam ";
	foreach my $chr(@{$sample{$sam}}){
	print "$chr ";	
	}	
	print "\n";
}
