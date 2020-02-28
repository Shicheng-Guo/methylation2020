#!/usr/bin/perl -w
#usage:
#calculateOnTarget.pl [BED] < pileup

use strict;

my $bed = $ARGV[0];
my %table;
my $window = 100;

my $countOnTarget=0;
my $target1x=0;
my $countAll=0;
my $targetSize=0;

open(BED, "$bed") || die("Error opening $bed\n");
while(my $line = <BED>){
	chomp($line);
	my @f = split /\t/, $line;
	my ($chr, $start, $end) = ($f[0], $f[1], $f[2]);
	$targetSize+=($end-$start+1);
	my $sbin = int($start/$window);
	my $ebin = int($end/$window);
	for(my $i = $sbin-2; $i<$ebin+2; $i++){
		push(@{$table{$i}}, "$chr-$start-$end");
	}
}
close(BED);

while(my $line = <STDIN>){
	chomp($line);
	my @f = split /\t/, $line;
	my ($chr, $pos, $cnt) = ($f[0], $f[1], $f[3]);
	$countAll+=$cnt;
	my $bin = int($pos/$window);
	next if(!$table{$bin});
	my @candidates = @{$table{$bin}};
	foreach my $can (@candidates){
		my ($qchr, $qstart, $qend) = split "-", $can;
		next if($qchr ne $chr);
		next if($pos > $qend || $pos < $qstart);
		$countOnTarget+=$cnt;
		$target1x++;
		last;
	}
}

print "On target bases = $countOnTarget\n";
print "Target size = $targetSize\n";
print "Specificity = $countOnTarget/$countAll = ", $countOnTarget/$countAll, "\n";
print "Sensitivity = $target1x/$targetSize = ", $target1x/$targetSize, "\n";
print "Enrichment = $countOnTarget/$target1x = ", $countOnTarget/$target1x, "\n";
print $countOnTarget/$countAll, "\t", $target1x/$targetSize, "\t", $countOnTarget/$target1x, "\n";



