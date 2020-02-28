#!/usr/bin/perl -w
#
use strict;

my %sampleFiles;
my %runInfo;
$runInfo{"MiSeq/PalssonLab/140630"}->{"dir"} = "140627_MISEQ_bam";
$runInfo{"Rapid run/IGM/140627"}->{"dir"} = "140627_SN1001_bam";
$runInfo{"Rapid run/IGM/140709"}->{"dir"} = "140709_SN1001_bam";
$runInfo{"Rapid run/ILMN/140712"}->{"dir"} = "140712_ILMN_bam";

$runInfo{"MiSeq/PalssonLab/140630"}->{"bam"} = "A0XX";
$runInfo{"Rapid run/IGM/140627"}->{"bam"} = "IndxXX";
$runInfo{"Rapid run/IGM/140709"}->{"bam"} = "IndxXX";
$runInfo{"Rapid run/ILMN/140712"}->{"bam"} = "IndxXX";

while(my $line = <STDIN>){
	chomp($line);
	next if($line =~ /TruSeq/);
	my @f = split /\t/, $line;
	my $run = $f[2];
	my $sampleId = $f[0];
	my $indexId = $f[1];
	my $laneNum = "s_".$f[3];
	my $files_dir = $runInfo{$run}->{"dir"};
	my $sampleBamId = $runInfo{$run}->{"bam"};
	$sampleBamId =~ s/XX/$indexId/g;
	#print $sampleBamId,"\n";
	my $files_list = `ls $files_dir`;
	my @files = split /\n/, $files_list;
	my $found = 0;
	foreach my $b (@files){
		chomp($b);
		if($run !~ m/MiSeq/ and $run !~ m/140709/){
			next if($b !~ m/$laneNum/);
		}
		my @parts = split /[-_\.]/, $b;
		foreach my $p (@parts){
			if($p eq $sampleBamId){
				push(@{$sampleFiles{$sampleId}}, "$files_dir/$b");
				$found = 1;
			}
		}
	}
	print $line, "\n" if($found == 0);
}

foreach my $sample_id (keys %sampleFiles){
	#print "#", $sample_id, "\n";
	my @files = @{$sampleFiles{$sample_id}};
	if(scalar(@files) == 1){
		print "cp ", $files[0], " $sample_id.merged.bam\n";
	}else{
		print "/home/dinh/softwares/samtools-0.1.18/samtools merge $sample_id.merged.bam ", join(" ", @files), "\n";
	}
}
