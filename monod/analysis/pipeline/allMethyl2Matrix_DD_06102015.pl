#!/usr/bin/perl  -w
#Should disable buffering of STDOUT
#Usage:
#allMethyl2Matrix.pl SAMPLES_LIST [minDepth] [minSamples] [sitesList/NA] [cnt/freq/depth]

#this script is more the new style methylFreq files

use strict;
use Statistics::Descriptive;

if(!$ARGV[0]){
	print "Usage:\n";
	print "allMethyl2Matrix.pl SAMPLES_LIST [minDepth] [minSamples] [sitesList/NA] [cnt/freq/depth]\n";
	print "\nthis script is more the new style methylFreq files, each site should have context either CG, CHG, or CHH\n";
	exit 0;
}


$| = 1;

my $list = $ARGV[0];
my $minSamples = 1;
$minSamples = $ARGV[2] if($ARGV[2]);
my $sitesList = "NA";
$sitesList = $ARGV[3] if($ARGV[3]);
my $style = "freq";
$style = $ARGV[4];

open(IN, "$list") || die("Error opening list file.");
my @fileList = <IN>;
close(IN);

my $minDepth=10;
$minDepth = $ARGV[1];
my @allSampleName;
my @allFiles;
my $print_header =1;

my %methylData;
my %methylCount;

if($sitesList eq "NA"){
	foreach my $file(@fileList){
		chomp($file);
		my @f = split /\t/, $file;
		my ($sampleID, $fileName) = ($f[0], $f[1]);
		open(INFILE, "$fileName") || die("Error in opening file $file.\n");
		while(my $line = <INFILE>){
			next if($line !~ /^chr/);
			chomp($line);
			my @fields = split(/\t/, $line);
			next if($fields[3] == 0);
			my $index = $fields[0] . ":" . $fields[1];
			$methylCount{$index}++;
		}
		close(INFILE);
		print "Finish scanning $file.\n";
	}
}else{
	print "Scanning $sitesList\n";
	open(INFILE, "$sitesList") || die("Error opening $sitesList\n");
	while(my $line = <INFILE>){
		chomp($line);
		$methylCount{$line}=$minSamples;
	}
	close(INFILE);
}

foreach my $file (@fileList){
	my @f = split /\t/, $file;
	my ($sampleID, $fileName) = ($f[0], $f[1]);
	open(INFILE, "$fileName") || die("Error in opening file $file.\n");	
	while(my $line = <INFILE>){
		next if($line !~ /^chr/);
		chomp($line);
		my @fields = split(/\t/, $line);
		my $index = $fields[0] . ":" . $fields[1];
		next if(!$methylCount{$index});
		next if($methylCount{$index} < $minSamples);
		next if($fields[3] == 0);
		my ($c, $t) = (0,0);
		for(my $i = 5; $i < scalar(@fields); $i+=2){
			$c = $fields[$i+1] if($fields[$i] eq 'C');
			$t = $fields[$i+1] if($fields[$i] eq 'T');
		}
		next if( ($c+$t)/$fields[3] < 0.9 );
		$methylData{$index}->{$sampleID}->{'C'}+=$c;
		$methylData{$index}->{$sampleID}->{'T'}+=$t;
	}
	close(INFILE);	
	print "Finish adding $file.\n";
	push(@allSampleName, $sampleID);
	push(@allFiles, $fileName);
}
undef(%methylCount);


my %fileCoverage;
my $outfile = "MethylMatrix.$list.$style";
open(OUT, ">$outfile") || die("Error writing $outfile.\n");
if($print_header){	
	print OUT "1chr_position\t", join("\t", @allSampleName), "\n";
}

foreach my $index(keys(%methylData)){			
		my @m;
		my $cnt =0;
		for(my $i=0; $i<scalar(@allSampleName); $i++){
			if($methylData{$index}->{$allSampleName[$i]}){
				my ($c,$t) = ($methylData{$index}->{$allSampleName[$i]}->{'C'}, $methylData{$index}->{$allSampleName[$i]}->{'T'});
				if($c+$t < $minDepth){
					push(@m, "NA");		
					next;
				}
				my $n = $c+$t;
				push(@m, sprintf("%4.3f", $c/($c+$t))) if($style eq "freq");
				push(@m, $c.":".$n) if($style eq "cnt");
				push(@m, $c+$t) if($style eq "depth");
				$cnt++;
				#$fileCoverage{$allFiles[$i]}->{"depth 1"}++ if($c + $t > 1);
				#$fileCoverage{$allFiles[$i]}->{"depth 5"}++ if($c + $t > 5);
				#$fileCoverage{$allFiles[$i]}->{"depth 10"}++ if($c + $t > 10);
			}else{
				push(@m, "NA");
			}
		}
		print OUT "$index\t", join("\t", @m), "\n" if($cnt >= $minSamples);
}	
close(OUT);

#foreach my $f(keys %fileCoverage){
#	print $f, "\t", $fileCoverage{$f}->{"depth 1"}, "\t", $fileCoverage{$f}->{"depth 5"}, "\t", $fileCoverage{$f}->{"depth 10"}, "\n";
#}
undef(%methylData);
undef(%fileCoverage);

